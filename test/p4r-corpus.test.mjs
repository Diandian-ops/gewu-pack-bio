/**
 * G0-RS-P4R-0 baseline corpus generation.
 *
 * The expensive test here runs real R. It is gated on the runtime being present,
 * and the gate is DECLARED rather than silent: an absent runtime produces a
 * visible skip with a reason, never a pass. A gate that reports success without
 * doing its work is the exact defect be0a6de fixed in the p2/p3/p4 acceptance
 * scripts, and it is not being reintroduced here.
 */

import { strict as assert } from 'node:assert';
import { execFileSync } from 'node:child_process';
import { existsSync, readFileSync, rmSync, statSync } from 'node:fs';
import { homedir } from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { fileURLToPath } from 'node:url';
import { findCoreRepo } from './helpers/core-repo.mjs';

import { fingerprintRealDataRoot } from '../../biof3-desktop/scripts/g0-rs-p4r-capture.mjs';
import { CORPUS_INPUTS, CORPUS_STEPS, generateCorpus } from '../../biof3-desktop/scripts/g0-rs-p4r-corpus.mjs';

const REPO_ROOT = findCoreRepo(); // BLOCK-83 P4：p4r 语料链跑在核心仓执行器上
const REPO_DB = path.join(REPO_ROOT, 'core', 'biof3-server', 'biof3.db');
const REAL_DATA_ROOT = path.join(homedir(), 'Documents', 'BioF3-Local');

/** Report why the R-backed test cannot run, or null when it can. */
function rRuntimeBlocker() {
  try {
    execFileSync('Rscript', ['-e', 'if(!requireNamespace("DESeq2",quietly=TRUE)) quit(status=1)'], {
      stdio: 'ignore',
      timeout: 60000,
    });
    return null;
  } catch {
    return 'requires Rscript with DESeq2 installed';
  }
}

const R_BLOCKER = rRuntimeBlocker();

test('the declared corpus inputs exist in the repository', () => {
  for (const [key, relative] of Object.entries(CORPUS_INPUTS)) {
    const full = path.join(REPO_ROOT, relative);
    assert.ok(existsSync(full), `${key} input missing: ${relative}`);
    assert.ok(statSync(full).size > 0, `${key} input is empty: ${relative}`);
  }
});

test('corpus steps are well formed and chain explicitly', () => {
  const ids = CORPUS_STEPS.map((step) => step.id);
  assert.equal(new Set(ids).size, ids.length, 'step ids must be unique');

  // heatmap is not a tool in this product; the plan's "heatmap" artifact comes
  // from deseq2's own output or from gene-correlation.
  assert.ok(!ids.includes('heatmap'), 'heatmap is not a real tool id');

  for (const step of CORPUS_STEPS) {
    assert.equal(typeof step.params, 'function', `${step.id} must build params programmatically`);
    if (step.dependsOn) {
      const upstreamIndex = CORPUS_STEPS.findIndex((s) => s.id === step.dependsOn);
      assert.ok(upstreamIndex !== -1, `${step.id} depends on unknown step ${step.dependsOn}`);
      assert.ok(
        upstreamIndex < CORPUS_STEPS.indexOf(step),
        `${step.id} must run after ${step.dependsOn}`,
      );
    }
  }

  // Nothing inspects a previous job, so a dependent step must name a concrete
  // upstream artifact path. Verify the wiring rather than trusting the comment.
  const volcano = CORPUS_STEPS.find((step) => step.id === 'volcano-plot');
  const params = volcano.params({
    repoRoot: REPO_ROOT,
    outputs: { deseq2: { 'deg_table.csv': '/somewhere/deg_table.csv' } },
  });
  assert.equal(params.deg_table, '/somewhere/deg_table.csv');

  // The enrichment observation is deliberately a real GO-BP run.  It must
  // not silently broaden into `all`, which would make this offline product
  // gate depend on KEGG availability and could turn a valid GO Artifact into
  // a failed Job merely because conditional KEGG outputs were absent.
  const enrichment = CORPUS_STEPS.find((step) => step.id === 'go-kegg');
  const enrichmentParams = enrichment.params({ repoRoot: REPO_ROOT, outputs: {} });
  assert.equal(enrichmentParams.analysis_type, 'GO-BP');
});

test('go-kegg declares conditional KEGG artifacts optional while keeping GO-BP truth required', () => {
  const definition = JSON.parse(readFileSync(
    path.join(REPO_ROOT, 'plugins/go-kegg/tool-definition.json'),
    'utf8',
  ));
  const outputs = new Map(definition.outputs.map((item) => [item.id, item]));
  assert.notEqual(outputs.get('go_bp_results')?.required, false);
  assert.notEqual(outputs.get('go_bp_dotplot')?.required, false);
  assert.equal(outputs.get('kegg_results')?.required, false);
  assert.equal(outputs.get('kegg_dotplot')?.required, false);
  assert.notEqual(outputs.get('report')?.required, false);
});

test(
  'generating the corpus runs the real chain and touches no real data',
  { timeout: 600000, skip: R_BLOCKER ?? false },
  async () => {
    const devDbBefore = existsSync(REPO_DB) ? statSync(REPO_DB).mtimeMs : null;
    const realRootBefore = fingerprintRealDataRoot(REAL_DATA_ROOT);

    const result = await generateCorpus({ keep: true, quiet: true });
    try {
      assert.equal(result.ok, true, `corpus failed: ${JSON.stringify(result.records)}`);
      assert.deepEqual(result.tables.missing, []);

      // Every step must have produced a real tool job and real artifacts.
      for (const record of result.records) {
        assert.equal(record.ok, true, `${record.toolId} did not succeed`);
        assert.match(record.toolJobId || '', /^tj_/, `${record.toolId} lacks a real tool job id`);
        assert.ok(record.artifactCount > 0, `${record.toolId} produced no artifacts`);
      }

      // The chain must be real: deseq2's DEG table and VST matrix are what the
      // downstream steps consumed, so they have to exist with real content.
      const deseq2 = result.outputs.deseq2 || {};
      for (const name of ['deg_table.csv', 'vst_matrix.csv', 'heatmap.png', 'report.html']) {
        const full = deseq2[name];
        assert.ok(full, `deseq2 did not emit ${name}`);
        assert.ok(statSync(full).size > 0, `${name} is empty`);
      }

      // Artifacts must live inside the disposable root, not anywhere real.
      for (const full of Object.values(deseq2)) {
        assert.ok(full.startsWith(result.root + path.sep), `artifact escaped the root: ${full}`);
      }
    } finally {
      rmSync(result.root, { recursive: true, force: true });
    }

    assert.equal(
      existsSync(REPO_DB) ? statSync(REPO_DB).mtimeMs : null,
      devDbBefore,
      'the developer database must be untouched by corpus generation',
    );
    assert.deepEqual(
      fingerprintRealDataRoot(REAL_DATA_ROOT),
      realRootBefore,
      'the real user data root must be untouched by corpus generation',
    );
  },
);

test('the R-backed corpus test is not silently skipped', () => {
  // If the runtime is missing the suite above skips, which node reports in its
  // summary. Assert the reason is stated so a skip can never read as coverage.
  if (R_BLOCKER) {
    assert.match(R_BLOCKER, /Rscript/, 'a skip must explain itself');
    console.error(`[p4r-corpus] real-chain coverage SKIPPED: ${R_BLOCKER}`);
  } else {
    assert.equal(R_BLOCKER, null);
  }
});
