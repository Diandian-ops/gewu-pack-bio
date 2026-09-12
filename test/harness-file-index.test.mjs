import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {
  createHarnessExecutionContext,
  createHarnessWorkflowPlan,
  createLocalHarnessPlatform,
  detectFile,
  diagnoseHarnessRunReadiness,
  indexHarnessScope,
  indexDirectory,
  inferColumnMappings,
  readIndexCache,
  rankHarnessToolsForRequest,
  repairHarnessParams,
  resolveHarnessScope,
  suggestInputBindings,
  summarizeIndexForPrompt,
  writeIndexCache,
} from '../scripts/harness-file-index.mjs';

function withTempDir(fn) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'biof3-harness-index-'));
  try {
    return fn(dir);
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
}

function readBuiltInManifest(id) {
  return JSON.parse(
    fs.readFileSync(path.join(process.cwd(), 'resources', 'built-in-plugins', id, 'tool-definition.json'), 'utf8'),
  );
}

test('detectFile detects DEG table headers and failure-repair useful hints', () => {
  withTempDir((dir) => {
    const file = path.join(dir, 'deg.csv');
    fs.writeFileSync(file, 'gene,log2FC,pvalue\nTP53,1.2,0.01\nEGFR,-2.1,0.001\n');
    const detected = detectFile(file);
    assert.equal(detected.kind, 'csv');
    assert.equal(detected.role, 'deg_table');
    assert.deepEqual(detected.table.headers, ['gene', 'log2FC', 'pvalue']);
    assert.equal(detected.table.hints.hasDegColumns, true);
    assert.equal(detected.table.columns[1].likelyType, 'number');
  });
});

test('detectFile infers rows_are_features expression matrix orientation', () => {
  withTempDir((dir) => {
    const file = path.join(dir, 'expr.csv');
    fs.writeFileSync(file, 'gene,S1,S2,S3\nA,1,2,3\nB,4,5,6\nC,7,8,9\n');
    const detected = detectFile(file);
    assert.equal(detected.role, 'expression_matrix');
    assert.equal(detected.table.matrixOrientation, 'rows_are_features');
    assert.equal(detected.table.columnCount, 4);
  });
});

test('detectFile infers rows_are_samples expression matrix orientation', () => {
  withTempDir((dir) => {
    const file = path.join(dir, 'sample_expr.csv');
    fs.writeFileSync(file, 'Sample,GENE001,GENE002,GENE003\nS1,1,2,3\nS2,4,5,6\nS3,7,8,9\n');
    const detected = detectFile(file);
    assert.equal(detected.role, 'expression_matrix');
    assert.equal(detected.table.matrixOrientation, 'rows_are_samples');
  });
});

test('indexDirectory builds metadata index and cache', () => {
  withTempDir((dir) => {
    fs.mkdirSync(path.join(dir, 'nested'));
    fs.writeFileSync(path.join(dir, 'groups.tsv'), 'sample\tgroup\nS1\tA\nS2\tB\n');
    fs.writeFileSync(path.join(dir, 'nested', 'notes.md'), '# Notes\nhello');
    const index = indexDirectory(dir);
    assert.equal(index.harnessIndexVersion, 1);
    assert.equal(index.fileCount, 2);
    assert.equal(index.byRole.sample_metadata, 1);
    assert.equal(index.byKind.text, 1);
    const cachePath = writeIndexCache(index, path.join(dir, '.biof3-harness-index.json'));
    const cached = readIndexCache(cachePath);
    assert.equal(cached.fileCount, 2);
    assert.ok(cached.entries.some((entry) => entry.relativePath === 'groups.tsv'));
  });
});

test('indexDirectory reuses unchanged cached entries and drops deleted files', () => {
  withTempDir((dir) => {
    fs.mkdirSync(path.join(dir, 'nested'));
    fs.writeFileSync(path.join(dir, 'keep.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    fs.writeFileSync(path.join(dir, 'nested', 'stale.tsv'), 'sample\tgroup\nS1\tA\n');

    const cachePath = path.join(dir, '.biof3-harness-index.json');
    const initial = indexDirectory(dir, { hash: false });
    writeIndexCache(initial, cachePath);

    const keepBefore = initial.entries.find((entry) => entry.relativePath === 'keep.csv');
    assert.ok(keepBefore, 'initial cache should contain keep.csv');

    fs.unlinkSync(path.join(dir, 'nested', 'stale.tsv'));
    fs.writeFileSync(path.join(dir, 'fresh.csv'), 'sample,group\nS1,A\n');

    const incremental = indexDirectory(dir, { hash: false, cache: cachePath });
    const keepAfter = incremental.entries.find((entry) => entry.relativePath === 'keep.csv');
    const freshAfter = incremental.entries.find((entry) => entry.relativePath === 'fresh.csv');

    assert.equal(incremental.cacheStats.hit, true);
    assert.equal(incremental.cacheStats.reusedEntries, 1);
    assert.equal(incremental.cacheStats.updatedEntries, 1);
    assert.equal(incremental.cacheStats.removedEntries, 1);
    assert.equal(incremental.fileCount, 2);
    assert.equal(keepAfter.inspectedAt, keepBefore.inspectedAt);
    assert.ok(freshAfter, 'fresh.csv should be indexed');
    assert.equal(freshAfter.relativePath, 'fresh.csv');
  });
});

test('suggestInputBindings ranks files for plugin inputs', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'demo_deg.csv'), 'gene,log2FC,pvalue\nTP53,1.2,0.01\n');
    fs.writeFileSync(path.join(dir, 'groups.tsv'), 'sample\tgroup\nS1\tA\n');
    const index = indexDirectory(dir, { hash: false });
    const suggestions = suggestInputBindings(
      {
        inputs: [
          { id: 'deg_table', label: '差异分析结果表', type: 'file', required: true, accept: '.csv' },
          { id: 'group_file', label: '样本分组文件', type: 'file', accept: '.csv,.tsv' },
        ],
      },
      index,
    );
    assert.equal(suggestions[0].selected.relativePath, 'demo_deg.csv');
    assert.equal(suggestions[1].selected.relativePath, 'groups.tsv');
    assert.ok(suggestions[0].selected.reasons.includes('role_deg_table'));
    assert.ok(suggestions[1].selected.reasons.includes('role_sample_metadata'));
  });
});

test('resolveHarnessScope resolves project path without database dependency', () => {
  withTempDir((dir) => {
    const resolved = resolveHarnessScope({ scope: 'project', projectId: 'p_demo', projectPath: dir });
    assert.equal(resolved.scope, 'project');
    assert.equal(resolved.id, 'p_demo');
    assert.equal(resolved.target, dir);
    assert.equal(resolved.platform, 'local');
    assert.equal(resolved.exists, true);
  });
});

test('createLocalHarnessPlatform resolves job and artifact paths behind adapter boundary', () => {
  withTempDir((dir) => {
    const jobsDir = path.join(dir, 'tool-jobs');
    const jobDir = path.join(jobsDir, 'job_adapter');
    fs.mkdirSync(path.join(jobDir, 'output'), { recursive: true });
    const platform = createLocalHarnessPlatform({ dataDir: dir, jobsDir });

    const jobScope = resolveHarnessScope({ platform, scope: 'job', jobId: 'job_adapter' });
    assert.equal(jobScope.platform, 'local');
    assert.equal(jobScope.partition, 'job');
    assert.equal(jobScope.target, jobDir);

    const artifactScope = resolveHarnessScope({ platform, scope: 'artifact', jobId: 'job_adapter' });
    assert.equal(artifactScope.partition, 'artifact');
    assert.equal(artifactScope.target, path.join(jobDir, 'output'));
  });
});

test('createLocalHarnessPlatform maps local paths to canonical execution paths', () => {
  withTempDir((dir) => {
    const projectDir = path.join(dir, 'project-a');
    fs.mkdirSync(projectDir, { recursive: true });
    fs.writeFileSync(path.join(projectDir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const platform = createLocalHarnessPlatform({
      dataDir: dir,
      pathMappings: [{ localRoot: dir, canonicalRoot: '/work/biof3' }],
    });

    const localFile = path.join(projectDir, 'deg.csv');
    assert.equal(platform.canonicalizePath(localFile), '/work/biof3/project-a/deg.csv');
    assert.equal(platform.resolveCanonicalPath('/work/biof3/project-a/deg.csv'), localFile);

    const scoped = indexHarnessScope({
      platform,
      scope: 'project',
      projectId: 'p_canonical',
      projectPath: projectDir,
      hash: false,
    });
    assert.equal(scoped.scope.canonicalTarget, '/work/biof3/project-a');
    assert.equal(scoped.index.canonicalRoot, '/work/biof3/project-a');
    assert.equal(scoped.index.entries[0].canonicalPath, '/work/biof3/project-a/deg.csv');
    assert.equal(scoped.promptSummary.canonicalRoot, '/work/biof3/project-a');
    assert.equal(scoped.promptSummary.entries[0].canonicalPath, '/work/biof3/project-a/deg.csv');
  });
});

test('createHarnessExecutionContext exposes runtime readiness and canonical workDir', () => {
  withTempDir((dir) => {
    const projectDir = path.join(dir, 'project-runtime');
    fs.mkdirSync(projectDir, { recursive: true });
    const platform = createLocalHarnessPlatform({
      dataDir: dir,
      pathMappings: [{ localRoot: dir, canonicalRoot: '/work/biof3' }],
      runtimeChecks: {
        r: { ready: true, command: '/opt/R/bin/Rscript' },
        python: false,
      },
    });

    const ctx = createHarnessExecutionContext({
      platform,
      scope: 'project',
      projectId: 'p_runtime',
      projectPath: projectDir,
      manifest: { language: 'r' },
    });

    assert.equal(ctx.platform, 'local');
    assert.equal(ctx.runtime, 'r');
    assert.equal(ctx.runtimeReady.ready, true);
    assert.equal(ctx.runtimeReady.command, '/opt/R/bin/Rscript');
    assert.equal(ctx.workDir, projectDir);
    assert.equal(ctx.canonicalWorkDir, '/work/biof3/project-runtime');
    assert.equal(platform.checkRuntimeReady('python').ready, false);
  });
});

test('createHarnessExecutionContext prefers manifest language over structured runtime metadata', () => {
  withTempDir((dir) => {
    const platform = createLocalHarnessPlatform({
      dataDir: dir,
      runtimeChecks: { r: { ready: true, command: '/opt/R/bin/Rscript' } },
    });
    const scopeInfo = platform.resolveScope({ scope: 'path', path: dir });
    const context = createHarnessExecutionContext({
      platform,
      scopeInfo,
      manifest: { language: 'r', runtime: { r: '>=4.3.0' } },
    });
    assert.equal(context.runtime, 'r');
    assert.equal(context.runtimeReady.command, '/opt/R/bin/Rscript');
  });
});

test('indexHarnessScope accepts custom platform adapter', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const platform = {
      id: 'test-platform',
      resolveScope() {
        return {
          scope: 'project',
          id: 'p_custom',
          target: dir,
          partition: 'project',
          platform: 'test-platform',
          exists: true,
        };
      },
    };
    const scoped = indexHarnessScope({ platform, scope: 'project', projectId: 'p_custom', hash: false });
    assert.equal(scoped.scope.platform, 'test-platform');
    assert.equal(scoped.index.byRole.deg_table, 1);
  });
});

test('indexHarnessScope indexes job and artifact partitions', () => {
  withTempDir((dir) => {
    const jobsDir = path.join(dir, 'tool-jobs');
    const jobDir = path.join(jobsDir, 'job_abc123');
    const outputDir = path.join(jobDir, 'output');
    fs.mkdirSync(outputDir, { recursive: true });
    fs.writeFileSync(path.join(jobDir, 'params.json'), '{"alpha":0.05}');
    fs.writeFileSync(path.join(outputDir, 'result.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');

    const jobIndexed = indexHarnessScope({
      scope: 'job',
      jobId: 'job_abc123',
      jobsDir,
      hash: false,
    });
    assert.equal(jobIndexed.scope.partition, 'job');
    assert.equal(jobIndexed.index.byRole.manifest || 0, 0);
    assert.ok(jobIndexed.index.entries.some((entry) => entry.relativePath === 'output/result.csv'));

    const artifactIndexed = indexHarnessScope({
      scope: 'artifact',
      jobId: 'job_abc123',
      jobsDir,
      hash: false,
    });
    assert.equal(artifactIndexed.scope.partition, 'artifact');
    assert.equal(artifactIndexed.index.fileCount, 1);
    assert.equal(artifactIndexed.index.entries[0].role, 'deg_table');
  });
});

test('summarizeIndexForPrompt removes bulky row samples', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'expr.csv'), 'gene,S1,S2,S3\nA,1,2,3\nB,4,5,6\n');
    const index = indexDirectory(dir, { hash: false });
    const summary = summarizeIndexForPrompt(index);
    assert.equal(summary.fileCount, 1);
    assert.equal(summary.entries[0].table.headers[0], 'gene');
    assert.equal(summary.entries[0].table.sampleRows, undefined);
    assert.equal(summary.entries[0].table.columns[0].name, 'gene');
  });
});

test('inferColumnMappings maps common DEG aliases', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg_alias.csv'), 'symbol,avg_log2FC,p_val\nTP53,1,0.01\n');
    const entry = detectFile(path.join(dir, 'deg_alias.csv'), { hash: false });
    const inferred = inferColumnMappings(entry, ['gene', 'log2FC', 'pvalue']);
    assert.deepEqual(inferred.mappings, {
      gene: 'symbol',
      log2FC: 'avg_log2FC',
      pvalue: 'p_val',
    });
    assert.deepEqual(inferred.missing, []);
  });
});

test('repairHarnessParams prefills file inputs and column mappings', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg_alias.csv'), 'symbol,avg_log2FC,p_val\nTP53,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const repaired = repairHarnessParams(
      {
        inputs: [
          { id: 'deg_table', label: '差异分析结果表', type: 'file', required: true, accept: '.csv' },
          { id: 'top_n', label: '标注 top N 基因', type: 'number', default: 10 },
        ],
      },
      index,
      { top_n: 999 },
    );
    assert.equal(repaired.repairedParams.deg_table, 'deg_alias.csv');
    assert.equal(repaired.repairedParams.top_n, 1);
    assert.deepEqual(repaired.repairedParams.columnMappings.deg_table, {
      gene: 'symbol',
      log2FC: 'avg_log2FC',
      pvalue: 'p_val',
    });
    assert.ok(repaired.repairs.some((repair) => repair.code === 'file_input_prefilled'));
    assert.ok(repaired.repairs.some((repair) => repair.code === 'number_param_clamped'));
  });
});

test('repairHarnessParams clamps PCA n_components to matrix capacity', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'expr.csv'), 'Sample,G1,G2\nS1,1,2\nS2,3,4\n');
    const index = indexDirectory(dir, { hash: false });
    const repaired = repairHarnessParams(
      {
        inputs: [
          { id: 'expr_matrix', label: '表达矩阵', type: 'file', required: true, accept: '.csv' },
          { id: 'n_components', label: '主成分数量', type: 'number', default: 2 },
        ],
      },
      index,
      { n_components: 10 },
    );
    assert.equal(repaired.repairedParams.expr_matrix, 'expr.csv');
    assert.equal(repaired.repairedParams.n_components, 2);
  });
});

test('repairHarnessParams preserves explicit column mappings and audits user value changes', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'symbol,avg_log2FC,p_val\nTP53,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const repaired = repairHarnessParams(volcano, index, {
      top_n: 999,
      columnMappings: { deg_table: { gene: 'symbol', log2FC: 'custom_fc' } },
    });

    assert.equal(repaired.repairedParams.columnMappings.deg_table.log2FC, 'custom_fc');
    assert.equal(repaired.repairedParams.columnMappings.deg_table.pvalue, 'p_val');
    assert.equal(repaired.audit.confirmationRequired, true);
    assert.deepEqual(repaired.audit.changedUserInputs, ['top_n']);
    const clamp = repaired.repairs.find((item) => item.code === 'number_param_clamped');
    assert.equal(clamp.origin, 'user_value_adjustment');
    assert.equal(clamp.confirmationRequired, true);
  });
});

test('weak-agent regression ranks volcano plot and repairs DEG params from index summary', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg_alias.csv'), 'symbol,avg_log2FC,p_val\nTP53,1,0.01\nEGFR,-2,0.001\n');
    fs.writeFileSync(path.join(dir, 'expr.csv'), 'Sample,G1,G2\nS1,1,2\nS2,3,4\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const pca = readBuiltInManifest('pca-explorer');
    const ranked = rankHarnessToolsForRequest('用当前项目里的差异分析结果画火山图，标注最多 999 个基因', [pca, volcano], index);

    assert.equal(ranked[0].toolId, 'volcano-plot');
    assert.ok(ranked[0].reasons.includes('intent_volcano'));
    assert.equal(
      ranked[0].inputSuggestions.find((binding) => binding.inputId === 'deg_table')?.selected?.relativePath,
      'deg_alias.csv',
    );

    const repaired = repairHarnessParams(volcano, index, { top_n: 999 });
    assert.equal(repaired.repairedParams.deg_table, 'deg_alias.csv');
    assert.equal(repaired.repairedParams.top_n, 2);
    assert.deepEqual(repaired.repairedParams.columnMappings.deg_table, {
      gene: 'symbol',
      log2FC: 'avg_log2FC',
      pvalue: 'p_val',
    });
  });
});

test('two-hop repair prioritizes the active upstream artifact and returns an executable absolute path', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'other_deg.csv'), 'gene,log2FC,pvalue\nA,1,0.01\n');
    fs.writeFileSync(path.join(dir, 'classified_table.csv'), 'gene,log2FC,pvalue\nB,2,0.001\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const inputSuggestions = suggestInputBindings(volcano, index, {
      preferredPaths: ['classified_table.csv'],
    });
    const selected = inputSuggestions.find((binding) => binding.inputId === 'deg_table')?.selected;
    assert.equal(selected.relativePath, 'classified_table.csv');
    assert.ok(selected.reasons.includes('active_artifact_preferred'));

    const repaired = repairHarnessParams(volcano, index, {}, {
      inputSuggestions,
      useAbsolutePaths: true,
    });
    assert.equal(repaired.repairedParams.deg_table, path.join(dir, 'classified_table.csv'));
  });
});

test('two-hop repair does not force an incompatible active artifact over a valid candidate', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'report.html'), '<html>report</html>');
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nA,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const inputSuggestions = suggestInputBindings(volcano, index, {
      preferredPaths: ['report.html'],
    });
    const selected = inputSuggestions.find((binding) => binding.inputId === 'deg_table')?.selected;

    assert.equal(selected.relativePath, 'deg.csv');
    assert.ok(!selected.reasons.includes('active_artifact_preferred'));
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, {}, {
      inputSuggestions,
      repair: repairHarnessParams(volcano, index, {}, { inputSuggestions, useAbsolutePaths: true }),
    });
    assert.equal(diagnostics.status, 'ready');
  });
});

test('weak-agent regression ranks PCA and binds expression matrix plus group file', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'expression_matrix.csv'), 'Sample,G1,G2,G3\nS1,1,2,3\nS2,4,5,6\nS3,7,8,9\n');
    fs.writeFileSync(path.join(dir, 'sample_groups.csv'), 'sample,group\nS1,A\nS2,B\nS3,A\n');
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const pca = readBuiltInManifest('pca-explorer');
    const ranked = rankHarnessToolsForRequest('对表达矩阵做 PCA 降维，并按分组文件上色', [volcano, pca], index);

    assert.equal(ranked[0].toolId, 'pca-explorer');
    assert.ok(ranked[0].reasons.includes('intent_pca'));
    assert.equal(
      ranked[0].inputSuggestions.find((binding) => binding.inputId === 'expr_matrix')?.selected?.relativePath,
      'expression_matrix.csv',
    );
    assert.equal(
      ranked[0].inputSuggestions.find((binding) => binding.inputId === 'group_file')?.selected?.relativePath,
      'sample_groups.csv',
    );

    const repaired = repairHarnessParams(pca, index, { n_components: 10 });
    assert.equal(repaired.repairedParams.expr_matrix, 'expression_matrix.csv');
    assert.equal(repaired.repairedParams.group_file, 'sample_groups.csv');
    assert.equal(repaired.repairedParams.n_components, 3);
  });
});

test('weak-agent regression ranks gene correlation heatmap and clamps top_n', () => {
  withTempDir((dir) => {
    fs.writeFileSync(
      path.join(dir, 'gene_expression.csv'),
      'Gene,S1,S2,S3,S4\nA,1,2,3,4\nB,2,3,4,5\nC,5,6,7,8\n',
    );
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const pca = readBuiltInManifest('pca-explorer');
    const correlation = readBuiltInManifest('gene-correlation');
    const ranked = rankHarnessToolsForRequest('用表达矩阵做基因相关性热图，top 999 个基因', [volcano, pca, correlation], index);

    assert.equal(ranked[0].toolId, 'gene-correlation');
    assert.ok(ranked[0].reasons.includes('intent_correlation'));
    assert.equal(
      ranked[0].inputSuggestions.find((binding) => binding.inputId === 'expr_matrix')?.selected?.relativePath,
      'gene_expression.csv',
    );

    const repaired = repairHarnessParams(correlation, index, { top_n: 999 });
    assert.equal(repaired.repairedParams.expr_matrix, 'gene_expression.csv');
    assert.equal(repaired.repairedParams.top_n, 3);
  });
});

test('diagnoseHarnessRunReadiness maps matrix orientation mismatch to failureMode issue', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'sample_rows.csv'), 'Sample,G1,G2,G3\nS1,1,2,3\nS2,4,5,6\nS3,7,8,9\n');
    const index = indexDirectory(dir, { hash: false });
    const correlation = readBuiltInManifest('gene-correlation');
    const diagnostics = diagnoseHarnessRunReadiness(correlation, index, {});

    assert.equal(diagnostics.status, 'needs_review');
    assert.equal(diagnostics.ok, true);
    const issue = diagnostics.issues.find((item) => item.code === 'matrix_orientation_mismatch');
    assert.ok(issue, '应识别 gene-correlation 的矩阵方向不匹配');
    assert.equal(issue.recoverable, true);
    assert.equal(issue.failureMode.code, 'matrix_orientation_mismatch');
    assert.equal(diagnostics.repairedParams.expr_matrix, 'sample_rows.csv');
  });
});

test('diagnoseHarnessRunReadiness blocks missing required file inputs', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'notes.md'), '# no usable data');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, {});

    assert.equal(diagnostics.ok, false);
    assert.equal(diagnostics.status, 'blocked');
    const issue = diagnostics.issues.find((item) => item.code === 'missing_required_file_inputs');
    assert.ok(issue, '应阻塞缺失 deg_table 的执行');
    assert.equal(issue.inputId, 'deg_table');
  });
});

test('diagnoseHarnessRunReadiness reports runtime_not_ready as recoverable blocker', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, {}, {
      runtimeReady: { runtime: 'r', ready: false, reason: 'r_runtime_not_found' },
    });

    assert.equal(diagnostics.ok, false);
    assert.equal(diagnostics.status, 'blocked');
    const issue = diagnostics.issues.find((item) => item.code === 'runtime_not_ready');
    assert.ok(issue, '应报告 runtime_not_ready');
    assert.equal(issue.recoverable, true);
    assert.match(issue.detail, /r_runtime_not_found/);
  });
});

test('diagnoseHarnessRunReadiness reports sample_group_mismatch for PCA inputs', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'expr.csv'), 'Sample,G1,G2,G3\nS1,1,2,3\nS2,4,5,6\nS3,7,8,9\n');
    fs.writeFileSync(path.join(dir, 'groups.csv'), 'sample,group\nX1,A\nX2,B\nX3,A\n');
    const index = indexDirectory(dir, { hash: false });
    const pca = readBuiltInManifest('pca-explorer');
    const diagnostics = diagnoseHarnessRunReadiness(pca, index, {});

    assert.equal(diagnostics.status, 'needs_review');
    const issue = diagnostics.issues.find((item) => item.code === 'sample_group_mismatch');
    assert.ok(issue, '应识别分组文件与表达矩阵不匹配');
    assert.equal(issue.recoverable, true);
    assert.equal(issue.severity, 'warning');
    assert.ok(diagnostics.nextActions.some((action) => action.id === 'align_group_samples'));
  });
});

test('diagnoseHarnessRunReadiness reports invalid_numeric_values and empty_significant_set for volcano inputs', () => {
  withTempDir((dir) => {
    fs.writeFileSync(
      path.join(dir, 'deg_bad.csv'),
      'gene,log2FC,pvalue\nTP53,foo,0.01\nEGFR,-2,bar\nBRCA1,0.3,0.9\n',
    );
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, { top_n: 10 });

    assert.equal(diagnostics.status, 'needs_review');
    assert.ok(diagnostics.issues.some((item) => item.code === 'invalid_numeric_values'));
  });
});

test('diagnostic gate blocks an unrepairable DEG schema and fully invalid numeric table', () => {
  withTempDir((dir) => {
    const volcano = readBuiltInManifest('volcano-plot');
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'name,effect,probability\nTP53,foo,bar\nEGFR,baz,qux\n');
    const index = indexDirectory(dir, { hash: false });
    const inputSuggestions = suggestInputBindings(volcano, index, { preferredPaths: ['deg.csv'] });
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, {}, { inputSuggestions });

    assert.equal(diagnostics.status, 'blocked');
    assert.equal(diagnostics.readiness, 'blocked');
    assert.ok(diagnostics.issues.some((item) => item.code === 'missing_required_columns'));
    assert.ok(diagnostics.issues.some((item) => item.code === 'invalid_numeric_values'));
  });
});

test('diagnoseHarnessRunReadiness reports empty_significant_set for volcano inputs', () => {
  withTempDir((dir) => {
    fs.writeFileSync(
      path.join(dir, 'deg_no_sig.csv'),
      'gene,log2FC,pvalue\nTP53,0.1,0.9\nEGFR,-0.2,0.8\nBRCA1,0.3,0.7\n',
    );
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const diagnostics = diagnoseHarnessRunReadiness(volcano, index, { top_n: 10, fc_threshold: 2, pval_threshold: 2 });

    assert.equal(diagnostics.status, 'needs_review');
    const empty = diagnostics.issues.find((item) => item.code === 'empty_significant_set');
    assert.ok(empty, '应识别没有满足阈值的显著基因');
    assert.equal(empty.recoverable, true);
  });
});

test('diagnoseHarnessRunReadiness surfaces next actions for blocked and ready states', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'notes.md'), '# no usable data');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const blocked = diagnoseHarnessRunReadiness(volcano, index, {});
    assert.equal(blocked.status, 'blocked');
    assert.ok(blocked.nextActions.some((action) => action.id === 'change_inputs'));
    assert.ok(blocked.nextActions.some((action) => action.id === 'inspect_summary'));

    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,1,0.01\n');
    const readyIndex = indexDirectory(dir, { hash: false });
    const ready = diagnoseHarnessRunReadiness(volcano, readyIndex, { top_n: 3 });
    assert.ok(ready.nextActions.some((action) => action.id === 'run'));
  });
});

test('workflow plan exposes deterministic ready lifecycle with canonical execution context', () => {
  withTempDir((dir) => {
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,2,0.001\n');
    const index = indexDirectory(dir, { hash: false });
    const volcano = readBuiltInManifest('volcano-plot');
    const platform = createLocalHarnessPlatform({
      localRoot: dir,
      canonicalRoot: '/workspace',
      runtimeChecks: { r: { ready: true, command: '/runtime/Rscript' } },
    });
    const scopeInfo = platform.resolveScope({ scope: 'path', path: dir });
    const executionContext = createHarnessExecutionContext({ platform, scopeInfo, manifest: volcano });
    const plan = createHarnessWorkflowPlan({
      requestText: '用差异结果画火山图',
      manifest: volcano,
      index,
      params: { top_n: 1 },
      executionContext,
      useAbsolutePaths: true,
    });

    assert.equal(plan.readiness, 'ready');
    assert.equal(plan.runnable, true);
    assert.equal(plan.confirmationRequired, false);
    assert.equal(plan.executionContext.canonicalWorkDir, '/workspace');
    assert.deepEqual(plan.stages.map((stage) => stage.id), [
      'inspect', 'infer', 'select', 'repair', 'confirm', 'run', 'verify', 'report',
    ]);
    assert.equal(plan.stages.find((stage) => stage.id === 'run').status, 'pending');
  });
});

test('workflow plan distinguishes recoverable review from blocked input', () => {
  withTempDir((dir) => {
    const volcano = readBuiltInManifest('volcano-plot');
    fs.writeFileSync(path.join(dir, 'deg.csv'), 'gene,log2FC,pvalue\nTP53,0.1,0.9\n');
    const recoverable = createHarnessWorkflowPlan({
      manifest: volcano,
      index: indexDirectory(dir, { hash: false }),
      params: { top_n: 1, fc_threshold: 2, pval_threshold: 2 },
    });
    assert.equal(recoverable.readiness, 'recoverable');
    assert.equal(recoverable.runnable, true);
    assert.equal(recoverable.confirmationRequired, true);
    assert.equal(recoverable.stages.find((stage) => stage.id === 'confirm').status, 'pending');

    fs.rmSync(path.join(dir, 'deg.csv'));
    fs.writeFileSync(path.join(dir, 'notes.md'), '# no table');
    const blocked = createHarnessWorkflowPlan({
      manifest: volcano,
      index: indexDirectory(dir, { hash: false }),
    });
    assert.equal(blocked.readiness, 'blocked');
    assert.equal(blocked.runnable, false);
    assert.equal(blocked.stages.find((stage) => stage.id === 'run').status, 'blocked');
  });
});
