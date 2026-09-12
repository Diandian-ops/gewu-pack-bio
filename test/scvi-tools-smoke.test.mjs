// test/scvi-tools-smoke.test.mjs — smoke tests for BLOCK-11 PHASE 1 SUB-1.5
//
// 验证（不依赖 scvi-tools/torch 真实运行）：
//   1. plugin.json / environment.json / acceptance.json 落盘且非空
//   2. plugin.json 通过 plugins/manifest.schema.json 校验
//   3. tool-definition.json 引用 8 输出文件
//   4. acceptance.json 四段 fixture (valid + missing-inputs + runtime-not-ready + empty-result)
//   5. placeholder PNG 存在且 sha256 在重写后稳定（worker C §0.5.5）
//   6. demo h5ad 落盘且能读 (counts layer + batch 字段)
//
// 注：本测试不调 ./scripts/scvi-tools.py 端到端（那是 dev Electron 边界内
// 单元测试范围——end-to-end 已由 SUB-1.5 worker 在 audit 阶段跑过）；
// 这里只验证 plugin 8 件套的合规与 schema 校验。
//
// 运行：
//   node --test test/scvi-tools-smoke.test.mjs

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
const pluginDir = resolve(repoRoot, 'plugins/scvi-tools');
const require = createRequire(import.meta.url);

test('plugin.json 存在 + 顶层 manifest.json 的 scvi-tools 条目通过 schema 校验', () => {
  const pluginPath = resolve(pluginDir, 'plugin.json');
  assert.ok(existsSync(pluginPath), 'plugin.json missing');
  assert.ok(statSync(pluginPath).size > 0, 'plugin.json empty');
  const plugin = JSON.parse(readFileSync(pluginPath, 'utf8'));
  assert.equal(plugin.id, 'scvi-tools');
  assert.equal(plugin.language, 'python');
  assert.ok(plugin.packages.python.scanpy, 'plugin.json must declare scanpy');
  assert.ok(plugin.packages.python.anndata, 'plugin.json must declare anndata');
  assert.ok(plugin.packages.python['scvi-tools'], 'plugin.json must declare scvi-tools');
  assert.ok(plugin.packages.python.torch, 'plugin.json must declare torch');
  assert.equal(plugin.device, 'cpu', 'plugin.json device must be cpu (G0 phase discipline)');

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
    manifest.plugins['scvi-tools'],
    'manifest.json scvi-tools entry missing after append'
  );
  const entry = manifest.plugins['scvi-tools'];
  // PluginEntry 必填: id/title/version/kind/artifacts
  assert.equal(entry.kind, 'plugin', 'kind must be plugin');
  assert.ok(entry.artifacts && typeof entry.artifacts === 'object', 'artifacts must be object');
  assert.ok(entry.artifacts['macos-arm64'], 'artifacts.macos-arm64 missing');
  assert.equal(entry.artifacts['macos-arm64'].url, 'built-in://scvi-tools');
  assert.ok(
    entry.packages.python['scvi-tools'],
    'manifest scvi-tools entry must declare scvi-tools package'
  );
  assert.equal(entry.device, 'cpu', 'manifest device must be cpu');
});

test('environment.json 存在且 packages.python 含 scanpy + anndata + scvi-tools + torch', () => {
  const path = resolve(pluginDir, 'environment.json');
  assert.ok(existsSync(path), 'environment.json missing');
  const env = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(env.packages && env.packages.python, 'environment.json packages.python missing');
  assert.ok(env.packages.python.scanpy, 'environment.json packages.python.scanpy missing');
  assert.ok(env.packages.python.anndata, 'environment.json packages.python.anndata missing');
  assert.ok(
    env.packages.python['scvi-tools'],
    'environment.json packages.python.scvi-tools missing'
  );
  assert.ok(env.packages.python.torch, 'environment.json packages.python.torch missing');
});

test('manifest.entry.json 存在且含 entryPoint 指向 scripts/scvi-tools.py', () => {
  const path = resolve(pluginDir, 'manifest.entry.json');
  assert.ok(existsSync(path), 'manifest.entry.json missing');
  const e = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(e.pluginId, 'scvi-tools', 'pluginId mismatch');
  assert.equal(
    e.entryPoint,
    'scripts/scvi-tools.py',
    'entryPoint should point to scripts/scvi-tools.py'
  );
  assert.ok(e.toolDefinitionFile, 'toolDefinitionFile missing');
  assert.ok(e.environmentFile, 'environmentFile missing');
  assert.equal(e.device, 'cpu', 'manifest.entry device must be cpu');
});

test('install.json 存在且声明 install-required + pip install 步骤', () => {
  const path = resolve(pluginDir, 'install.json');
  assert.ok(existsSync(path), 'install.json missing');
  const i = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(i.schemaVersion, 1, 'install.json schemaVersion mismatch');
  assert.equal(i.pluginId, 'scvi-tools');
  assert.equal(
    i.policy,
    'install-required',
    'install.json policy must be install-required (scvi-tools/torch not in biof3-py-runtime)'
  );
  assert.ok(Array.isArray(i.steps), 'install.steps should be array');
  assert.ok(i.steps.length > 0, 'install.steps must have at least 1 pip install step');
  const step = i.steps[0];
  assert.equal(step.kind, 'pip-install', 'install step kind must be pip-install');
  assert.ok(
    step.command.includes('scvi-tools') && step.command.includes('torch'),
    'install step command must install both scvi-tools and torch'
  );
});

test('scripts/scvi-tools.py 存在且为合法 Python 入口', () => {
  const path = resolve(pluginDir, 'scripts/scvi-tools.py');
  assert.ok(existsSync(path), 'scripts/scvi-tools.py missing');
  const stats = statSync(path);
  assert.ok(stats.size > 1000, `scripts/scvi-tools.py too small (${stats.size} bytes)`);
  const body = readFileSync(path, 'utf8');
  assert.ok(body.startsWith('#!/usr/bin/env python3'), 'should have shebang');
  assert.ok(
    body.includes('sys.argv[1]') || body.includes('commandArgs'),
    'should honor BioF3 dispatcher contract'
  );
  assert.ok(
    body.includes('_canonicalize_png') && body.includes('Pillow'),
    'should canonicalize PNG via Pillow (worker C §0.5.5)'
  );
  assert.ok(
    body.includes('scvi.model.SCVI'),
    'should call scvi.model.SCVI for deep generative model'
  );
  assert.ok(
    body.includes('get_latent_representation'),
    'should extract latent via get_latent_representation'
  );
  assert.ok(
    body.includes('get_normalized_expression'),
    'should extract normalized expression'
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
  assert.ok(Array.isArray(input.required) && input.required.includes('anndata_file'));
  assert.ok(Array.isArray(output.required) && output.required.includes('outputs'));
});

test('tool-definition.json 含 8 个 outputs 且 ipcChannel + failureModes 完整', () => {
  const path = resolve(pluginDir, 'tool-definition.json');
  assert.ok(existsSync(path), 'tool-definition.json missing');
  const td = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(td.id, 'scvi-tools');
  assert.equal(td.language, 'python');
  assert.equal(td.category, 'single-cell');
  assert.ok(
    Array.isArray(td.outputs) && td.outputs.length === 8,
    `tool-definition.json outputs length must be 8, got ${td.outputs?.length}`
  );

  const expectedOutputs = [
    'latent_embedding',
    'normalized_expression',
    'integrated_h5ad',
    'latent_umap',
    'latent_umap_batch',
    'training_loss',
    'report',
    'manifest',
  ];
  const ids = td.outputs.map(o => o.id);
  for (const id of expectedOutputs) {
    assert.ok(ids.includes(id), `outputs missing: ${id}`);
  }

  assert.ok(td.ipcChannel === 'tools:invoke/scvi-tools', 'ipcChannel mismatch');
  assert.ok(
    td.failureModes && td.failureModes.length === 3,
    `failureModes should have 3 entries, got ${td.failureModes?.length}`
  );
  const failureCodes = td.failureModes.map(f => f.code);
  assert.ok(failureCodes.includes('missing_required_file_inputs'));
  assert.ok(failureCodes.includes('runtime_not_ready'));
  assert.ok(failureCodes.includes('empty_result'));
});

test('acceptance.json 4 段 fixture (valid + missing-inputs + runtime-not-ready + empty-result)', () => {
  const path = resolve(pluginDir, 'acceptance.json');
  assert.ok(existsSync(path), 'acceptance.json missing');
  const a = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(a.pluginId, 'scvi-tools');
  assert.ok(Array.isArray(a.fixtures));
  const ids = a.fixtures.map(f => f.id);
  assert.ok(ids.includes('valid'), 'fixture valid missing');
  assert.ok(ids.includes('missing-inputs'), 'fixture missing-inputs missing');
  assert.ok(ids.includes('runtime-not-ready'), 'fixture runtime-not-ready missing');
  assert.ok(ids.includes('empty-result'), 'fixture empty-result missing');

  // runtime-not-ready 必须显式声明 skipOnSupportedPlatform: true
  const rnr = a.fixtures.find(f => f.id === 'runtime-not-ready');
  assert.equal(
    rnr.skipOnSupportedPlatform,
    true,
    'runtime-not-ready fixture must set skipOnSupportedPlatform: true'
  );

  // runtimeExtension 必须声明 scvi-tools + torch 额外依赖
  assert.ok(a.runtimeExtension, 'runtimeExtension must be declared');
  assert.ok(
    Array.isArray(a.runtimeExtension.extraDependencies) &&
      a.runtimeExtension.extraDependencies.length >= 2,
    'runtimeExtension must declare scvi-tools + torch extraDependencies'
  );
  const depNames = a.runtimeExtension.extraDependencies.map(d => d.name);
  assert.ok(depNames.includes('scvi-tools'), 'extraDependencies must include scvi-tools');
  assert.ok(depNames.includes('torch'), 'extraDependencies must include torch');
});

test('placeholder PNG 存在且 > 200 bytes（hash rewrite 后稳定输出）', () => {
  const path = resolve(pluginDir, 'tools/figures/latent_umap_placeholder.png');
  assert.ok(existsSync(path), 'tools/figures/latent_umap_placeholder.png missing');
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

test('demo h5ad 存在且 > 100KB', () => {
  const path = resolve(pluginDir, 'demo-data/pbmc_3k_mini.h5ad');
  assert.ok(existsSync(path), 'demo-data/pbmc_3k_mini.h5ad missing');
  const stats = statSync(path);
  assert.ok(stats.size > 100000, `demo h5ad too small: ${stats.size} bytes`);
});

test('顶层 manifest.json 已包含 scvi-tools 条目', () => {
  const path = resolve(repoRoot, 'plugins/manifest.json');
  const manifest = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(manifest.plugins['scvi-tools'], 'manifest.json scvi-tools entry missing');
  assert.equal(manifest.plugins['scvi-tools'].id, 'scvi-tools');
  assert.equal(manifest.plugins['scvi-tools'].language, 'python');
  assert.equal(manifest.plugins['scvi-tools'].kind, 'plugin');
});

test('scvi-tools plugin 目录结构完整（8 件套文件全部存在）', () => {
  const required = [
    'plugin.json',
    'environment.json',
    'tool-definition.json',
    'acceptance.json',
    'manifest.entry.json',
    'install.json',
    'README.md',
    'scvi-tools-report-template.html',
    'scripts/scvi-tools.py',
    'schema/input.schema.json',
    'schema/output.schema.json',
    'demo-data/pbmc_3k_mini.h5ad',
    'demo-data/build_demo_h5ad.py',
    'tools/figures/latent_umap_placeholder.png',
    'tools/figures/generate_placeholder_figure.py',
  ];
  const missing = required.filter(rel => !existsSync(resolve(pluginDir, rel)));
  assert.equal(missing.length, 0, `missing files: ${missing.join(', ')}`);
});

test('scvi-tools plugin 目录无未提交孤儿文件', () => {
  const requiredDirs = ['scripts', 'schema', 'demo-data', 'tools/figures', 'tools'];
  for (const d of requiredDirs) {
    const dir = resolve(pluginDir, d);
    assert.ok(existsSync(dir), `directory missing: ${d}`);
    const items = readdirSync(dir);
    assert.ok(items.length > 0, `directory empty: ${d}`);
  }
});
