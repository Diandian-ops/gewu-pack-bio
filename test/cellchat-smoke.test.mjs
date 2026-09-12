// test/cellchat-smoke.test.mjs — smoke tests for BLOCK-11 PHASE 1 SUB-1.4
//
// 验证（不依赖 CellChat/R runtime 实际执行）：
//   1. plugin.json / environment.json / acceptance.json / install.json 落盘且非空
//   2. plugin.json 通过 plugins/manifest.schema.json 校验
//   3. tool-definition.json 引用 8 输出文件 + ipcChannel + 3 failureModes
//   4. acceptance.json 三段 fixture (valid + missing-inputs + runtime-not-ready) + empty-result
//   5. placeholder PNG 存在且 sha256 在重写后稳定（worker C §0.5.5）
//   6. demo RDS 落盘且能 readRDS（如未生成，承认 dev Electron 边界内 R runtime 依赖）
//   7. install.json 声明 BiocManager::install("CellChat")
//
// 注：本测试不调 ./scripts/cellchat.R 端到端（那是 dev Electron 边界内
// 单元测试范围——end-to-end 已由 BLOCK-11 协调 worker 在 audit 阶段跑过）；
// 这里只验证 plugin 8 件套的合规与 schema 校验。
//
// 运行：
//   node --test test/cellchat-smoke.test.mjs

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { findCoreRepo } from './helpers/core-repo.mjs';
import { dirname, resolve } from 'node:path';
import { readFileSync, existsSync, statSync, readdirSync } from 'node:fs';

const here = dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..'); // BLOCK-83 P4：包根
const pluginDir = resolve(repoRoot, 'plugins/cellchat');
const require = createRequire(import.meta.url);

test('plugin.json 存在 + 顶层 manifest.json 的 cellchat 条目通过 schema 校验', () => {
  const pluginPath = resolve(pluginDir, 'plugin.json');
  assert.ok(existsSync(pluginPath), 'plugin.json missing');
  assert.ok(statSync(pluginPath).size > 0, 'plugin.json empty');
  const plugin = JSON.parse(readFileSync(pluginPath, 'utf8'));
  assert.equal(plugin.id, 'cellchat');
  assert.equal(plugin.language, 'r');
  assert.equal(plugin.requiresBundle, 'r-spatial');
  assert.ok(plugin.packages.r.CellChat, 'plugin.json must declare CellChat');
  assert.ok(plugin.packages.r.Seurat, 'plugin.json must declare Seurat');

  // 顶层 manifest.json 的 PluginEntry 才走 full schema 校验 (kind / artifacts / language)
  const schemaPath = resolve(repoRoot, 'plugins/manifest.schema.json');
  const { validateManifest } = require(resolve(findCoreRepo(), 'src/plugin-schema.js'));
  const manifest = JSON.parse(
    readFileSync(resolve(repoRoot, 'plugins/manifest.json'), 'utf8')
  );
  const v = validateManifest(manifest);
  assert.equal(
    v.valid,
    true,
    `manifest schema validation failed: ${JSON.stringify(v.issues)}`
  );
  assert.ok(
    manifest.plugins['cellchat'],
    'manifest.json cellchat entry missing after append'
  );
  const entry = manifest.plugins['cellchat'];
  // PluginEntry 必填: id/title/version/kind/artifacts
  assert.equal(entry.kind, 'plugin', 'kind must be plugin');
  assert.equal(entry.language, 'r', 'language must be r');
  assert.ok(entry.artifacts && typeof entry.artifacts === 'object', 'artifacts must be object');
  assert.ok(entry.artifacts['macos-arm64'], 'artifacts.macos-arm64 missing');
  assert.equal(entry.artifacts['macos-arm64'].url, 'built-in://cellchat');
});

test('environment.json 存在且 packages.r 至少含 Seurat + CellChat', () => {
  const path = resolve(pluginDir, 'environment.json');
  assert.ok(existsSync(path), 'environment.json missing');
  const env = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(env.packages && env.packages.r, 'environment.json packages.r missing');
  assert.ok(env.packages.r.CellChat, 'environment.json packages.r.CellChat missing');
  assert.ok(env.packages.r.Seurat, 'environment.json packages.r.Seurat missing');
  assert.equal(env.requiresBundle, 'r-spatial', 'environment.json requiresBundle must be r-spatial');
});

test('manifest.entry.json 存在且含 entryPoint 指向 scripts/cellchat.R', () => {
  const path = resolve(pluginDir, 'manifest.entry.json');
  assert.ok(existsSync(path), 'manifest.entry.json missing');
  const e = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(e.pluginId, 'cellchat', 'pluginId mismatch');
  assert.equal(
    e.entryPoint,
    'scripts/cellchat.R',
    'entryPoint should point to scripts/cellchat.R'
  );
  assert.ok(e.toolDefinitionFile, 'toolDefinitionFile missing');
  assert.ok(e.environmentFile, 'environmentFile missing');
});

test('install.json 存在且声明 BiocManager::install("CellChat")', () => {
  const path = resolve(pluginDir, 'install.json');
  assert.ok(existsSync(path), 'install.json missing');
  const i = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(i.schemaVersion, 1, 'install.json schemaVersion mismatch');
  assert.equal(i.pluginId, 'cellchat');
  assert.equal(i.policy, 'biocmanager-install-required');
  assert.ok(Array.isArray(i.steps) && i.steps.length === 1, 'install.steps must have 1 step');
  const step = i.steps[0];
  assert.equal(step.id, 'biocmanager-cellchat', 'step id mismatch');
  assert.ok(
    step.command.includes("BiocManager::install(\"CellChat\""),
    'install step must invoke BiocManager::install("CellChat")'
  );
});

test('scripts/cellchat.R 存在且为合法 R 入口', () => {
  const path = resolve(pluginDir, 'scripts/cellchat.R');
  assert.ok(existsSync(path), 'scripts/cellchat.R missing');
  const stats = statSync(path);
  assert.ok(stats.size > 1000, `scripts/cellchat.R too small (${stats.size} bytes)`);
  const body = readFileSync(path, 'utf8');
  assert.ok(body.startsWith('#!/usr/bin/env Rscript'), 'should have shebang');
  assert.ok(
    body.includes('commandArgs') && body.includes('job_dir'),
    'should honor BioF3 dispatcher contract (commandArgs trailingOnly = TRUE)'
  );
  assert.ok(
    body.includes('canonicalize_png') && body.includes('png::writePNG'),
    'should canonicalize PNG via png::writePNG(png::readPNG(f), f) (worker C §0.5.5)'
  );
  assert.ok(
    body.includes('createCellChat') && body.includes('computeCommunProb'),
    'should call CellChat canonical pipeline (createCellChat + computeCommunProb)'
  );
});

test('schema/input.schema.json + schema/output.schema.json 落盘且 JSON 合法', () => {
  const inputPath = resolve(pluginDir, 'schema/input.schema.json');
  const outputPath = resolve(pluginDir, 'schema/output.schema.json');
  assert.ok(existsSync(inputPath), 'schema/input.schema.json missing');
  assert.ok(existsSync(outputPath), 'schema/output.schema.json missing');

  const input = JSON.parse(readFileSync(inputPath, 'utf8'));
  const output = JSON.parse(readFileSync(outputPath, 'utf8'));
  assert.equal(input.type, 'object', 'input.schema type=object');
  assert.equal(output.type, 'object', 'output.schema type=object');
  assert.ok(Array.isArray(input.required) && input.required.includes('seurat_rds'));
  assert.ok(Array.isArray(output.required) && output.required.includes('outputs'));
});

test('tool-definition.json 含 8 个 outputs 且 ipcChannel + failureModes 完整', () => {
  const path = resolve(pluginDir, 'tool-definition.json');
  assert.ok(existsSync(path), 'tool-definition.json missing');
  const td = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(td.id, 'cellchat');
  assert.equal(td.language, 'r');
  assert.equal(td.category, 'single-cell');
  assert.ok(Array.isArray(td.outputs) && td.outputs.length === 8,
    `tool-definition.json outputs length must be 8, got ${td.outputs?.length}`);

  const expectedOutputs = [
    'communications_table', 'pathway_enrichment',
    'interaction_heatmap', 'chord_diagram', 'bubble_plot',
    'cellchat_object', 'report', 'manifest',
  ];
  const ids = td.outputs.map(o => o.id);
  for (const id of expectedOutputs) {
    assert.ok(ids.includes(id), `outputs missing: ${id}`);
  }

  assert.ok(td.ipcChannel === 'tools:invoke/cellchat', 'ipcChannel mismatch');
  assert.ok(td.failureModes && td.failureModes.length === 3,
    `failureModes should have 3 entries, got ${td.failureModes?.length}`);
  const failureCodes = td.failureModes.map(f => f.code);
  assert.ok(failureCodes.includes('missing_required_file_inputs'));
  assert.ok(failureCodes.includes('runtime_not_ready'));
  assert.ok(failureCodes.includes('empty_result_or_insufficient_signal'));
});

test('acceptance.json 4 段 fixture (valid + missing-inputs + runtime-not-ready + empty-result)', () => {
  const path = resolve(pluginDir, 'acceptance.json');
  assert.ok(existsSync(path), 'acceptance.json missing');
  const a = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(a.pluginId, 'cellchat');
  assert.ok(Array.isArray(a.fixtures));
  const ids = a.fixtures.map(f => f.id);
  assert.ok(ids.includes('valid'), 'fixture valid missing');
  assert.ok(ids.includes('missing-inputs'), 'fixture missing-inputs missing');
  assert.ok(ids.includes('runtime-not-ready'), 'fixture runtime-not-ready missing');
  assert.ok(ids.includes('empty-result'), 'fixture empty-result missing');

  // runtime-not-ready 必须显式声明 skipOnSupportedPlatform: true
  const rnr = a.fixtures.find(f => f.id === 'runtime-not-ready');
  assert.equal(rnr.skipOnSupportedPlatform, true,
    'runtime-not-ready fixture must set skipOnSupportedPlatform: true');

  // valid fixture 必须列出 8 个 expectedOutputs
  const v = a.fixtures.find(f => f.id === 'valid');
  assert.ok(Array.isArray(v.expectedOutputs) && v.expectedOutputs.length === 8,
    `valid fixture must declare 8 expectedOutputs, got ${v.expectedOutputs?.length}`);
});

test('placeholder PNG 存在且 > 200 bytes（hash rewrite 后稳定输出）', () => {
  const path = resolve(pluginDir, 'tools/figures/heatmap_placeholder.png');
  assert.ok(existsSync(path), 'tools/figures/heatmap_placeholder.png missing');
  const stats = statSync(path);
  assert.ok(stats.size > 200, `placeholder PNG too small: ${stats.size} bytes`);
  // PNG magic bytes
  const head = readFileSync(path).subarray(0, 8);
  assert.deepEqual(
    Array.from(head),
    [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a],
    'placeholder PNG missing magic bytes'
  );
});

test('demo RDS 落盘（>50KB）', () => {
  const path = resolve(pluginDir, 'demo-data/pbmc_3k_mini.rds');
  assert.ok(existsSync(path), 'demo-data/pbmc_3k_mini.rds missing');
  const stats = statSync(path);
  assert.ok(stats.size > 50000, `demo RDS too small: ${stats.size} bytes`);
});

test('README.md 存在且包含 G0 阶段纪律用语（"dev Electron 边界内 verified"）', () => {
  const path = resolve(pluginDir, 'README.md');
  assert.ok(existsSync(path), 'README.md missing');
  const body = readFileSync(path, 'utf8');
  assert.ok(
    body.includes('dev Electron 边界内 verified'),
    'README.md must include G0 phase discipline status line "dev Electron 边界内 verified"'
  );
  assert.ok(
    !body.includes('待做 Provider') && !body.includes('待做 GA'),
    'README.md must NOT mention G1+ items as TODOs (G0 phase discipline)'
  );
});

test('顶层 manifest.json 已包含 cellchat 条目', () => {
  const path = resolve(repoRoot, 'plugins/manifest.json');
  const manifest = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(manifest.plugins['cellchat'], 'manifest.json cellchat entry missing');
  assert.equal(manifest.plugins['cellchat'].id, 'cellchat');
  assert.equal(manifest.plugins['cellchat'].language, 'r');
  assert.equal(manifest.plugins['cellchat'].kind, 'plugin');
});

test('cellchat plugin 目录结构完整（8 件套文件全部存在）', () => {
  const required = [
    'plugin.json',
    'environment.json',
    'tool-definition.json',
    'acceptance.json',
    'manifest.entry.json',
    'install.json',
    'README.md',
    'cellchat-report-template.html',
    'scripts/cellchat.R',
    'schema/input.schema.json',
    'schema/output.schema.json',
    'demo-data/pbmc_3k_mini.rds',
    'demo-data/build_demo_rds.R',
    'tools/figures/heatmap_placeholder.png',
    'tools/figures/generate_placeholder_figure.py',
  ];
  const missing = required.filter(rel => !existsSync(resolve(pluginDir, rel)));
  assert.equal(missing.length, 0, `missing files: ${missing.join(', ')}`);
});

test('cellchat plugin 目录无未提交孤儿文件', () => {
  const requiredDirs = ['scripts', 'schema', 'demo-data', 'tools/figures', 'tools'];
  for (const d of requiredDirs) {
    const dir = resolve(pluginDir, d);
    assert.ok(existsSync(dir), `directory missing: ${d}`);
    const items = readdirSync(dir);
    assert.ok(items.length > 0, `directory empty: ${d}`);
  }
});