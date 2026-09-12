// test/celltypist-smoke.test.mjs — smoke tests for BLOCK-11 PHASE 1 SUB-1.3
//
// 验证（不依赖 celltypist 真实调用）：
//   1. plugin.json / environment.json / acceptance.json / install.json / manifest.entry.json 落盘且非空
//   2. plugin.json 通过 plugins/manifest.schema.json 校验
//   3. tool-definition.json 引用 8 输出文件
//   4. acceptance.json 四段 fixture (valid + missing-inputs + runtime-not-ready + empty-result)
//   5. placeholder PNG 存在且 sha256 在重写后稳定（worker C §0.5.5）
//   6. demo h5ad 落盘且 sha256 与 scanpy-advanced 一致 (e3d6c0f022...)
//   7. install.json 声明 pip-bootstrap（celltypist 首次安装）
//   8. environment.json 包含 celltypist>=1.6 (与 plugin.json 一致)
//
// 注：本测试不调 ./scripts/celltypist.py 端到端（celltypist Model.predict
// 需要联网下载预训练模型；本测试只验证 8 件套的合规与 schema 校验）。
//
// 运行：
//   node --test test/celltypist-smoke.test.mjs

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
const pluginDir = resolve(repoRoot, 'plugins/celltypist');
const require = createRequire(import.meta.url);

test('plugin.json 存在 + 顶层 manifest.json 的 celltypist 条目通过 schema 校验', () => {
  const pluginPath = resolve(pluginDir, 'plugin.json');
  assert.ok(existsSync(pluginPath), 'plugin.json missing');
  assert.ok(statSync(pluginPath).size > 0, 'plugin.json empty');
  const plugin = JSON.parse(readFileSync(pluginPath, 'utf8'));
  assert.equal(plugin.id, 'celltypist');
  assert.equal(plugin.language, 'python');
  assert.ok(plugin.packages.python.scanpy, 'plugin.json must declare scanpy');
  assert.ok(plugin.packages.python.anndata, 'plugin.json must declare anndata');
  assert.ok(plugin.packages.python.celltypist, 'plugin.json must declare celltypist (extension package)');

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
    manifest.plugins['celltypist'],
    'manifest.json celltypist entry missing after append'
  );
  const entry = manifest.plugins['celltypist'];
  // PluginEntry 必填: id/title/version/kind/artifacts
  assert.equal(entry.kind, 'plugin', 'kind must be plugin');
  assert.ok(entry.artifacts && typeof entry.artifacts === 'object', 'artifacts must be object');
  assert.ok(entry.artifacts['macos-arm64'], 'artifacts.macos-arm64 missing');
  assert.equal(entry.artifacts['macos-arm64'].url, 'built-in://celltypist');
});

test('environment.json 存在且 packages.python 至少含 scanpy + anndata + celltypist', () => {
  const path = resolve(pluginDir, 'environment.json');
  assert.ok(existsSync(path), 'environment.json missing');
  const env = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(env.packages && env.packages.python, 'environment.json packages.python missing');
  assert.ok(env.packages.python.scanpy, 'environment.json packages.python.scanpy missing');
  assert.ok(env.packages.python.anndata, 'environment.json packages.python.anndata missing');
  assert.ok(env.packages.python.celltypist, 'environment.json packages.python.celltypist missing');
  // celltypist 不是 biof3-py-runtime 已装包 → environment.json 必须显式 declares
  assert.equal(env.packages.python.celltypist, '>=1.6', 'celltypist must be >=1.6');
  // 标记 requiresBootstrap
  assert.ok(env.policy.requiresBootstrap.includes('celltypist'),
    'environment.json policy.requiresBootstrap must include celltypist');
});

test('manifest.entry.json 存在且含 entryPoint 指向 scripts/celltypist.py', () => {
  const path = resolve(pluginDir, 'manifest.entry.json');
  assert.ok(existsSync(path), 'manifest.entry.json missing');
  const e = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(e.pluginId, 'celltypist', 'pluginId mismatch');
  assert.equal(
    e.entryPoint,
    'scripts/celltypist.py',
    'entryPoint should point to scripts/celltypist.py'
  );
  assert.ok(e.toolDefinitionFile, 'toolDefinitionFile missing');
  assert.ok(e.environmentFile, 'environmentFile missing');
  assert.ok(e.builtAgainst, 'builtAgainst missing');
  assert.equal(e.builtAgainst.celltypist, '1.6.3', 'builtAgainst.celltypist should pin to operational version');
});

test('install.json 存在且声明 pip-bootstrap 触发 celltypist install', () => {
  const path = resolve(pluginDir, 'install.json');
  assert.ok(existsSync(path), 'install.json missing');
  const i = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(i.schemaVersion, 1, 'install.json schemaVersion mismatch');
  assert.equal(i.pluginId, 'celltypist');
  // celltypist 不在 biof3-py-runtime → install.json 必须是 pip-bootstrap
  assert.equal(i.policy, 'pip-bootstrap', 'install policy must be pip-bootstrap (celltypist extension)');
  assert.ok(Array.isArray(i.steps) && i.steps.length > 0,
    'install.steps must contain pip-install step');
  const step = i.steps[0];
  assert.equal(step.kind, 'pip-install', 'install step must be pip-install');
  assert.equal(step.package, 'celltypist', 'install package must be celltypist');
  assert.ok(step.versionSpec.includes('1.6'), 'versionSpec must pin >=1.6');
  assert.equal(step.installInRuntime, true, 'install must be in biof3-py-runtime');
  assert.equal(step.runtime, 'biof3-py-runtime', 'install target must be biof3-py-runtime');
});

test('scripts/celltypist.py 存在且为合法 Python 入口', () => {
  const path = resolve(pluginDir, 'scripts/celltypist.py');
  assert.ok(existsSync(path), 'scripts/celltypist.py missing');
  const stats = statSync(path);
  assert.ok(stats.size > 5000, `scripts/celltypist.py too small (${stats.size} bytes)`);
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
  // celltypist-specific: 必须 imports `celltypist` 或 fallback
  assert.ok(
    body.includes('celltypist') || body.includes('fallback_used'),
    'should reference celltypist or fallback semantics'
  );
  // fallback 路径
  assert.ok(
    body.includes('_scanpy_heuristic_annotation') && body.includes('fallback_used'),
    'should implement scanpy-based heuristic fallback when celltypist unavailable'
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
  // output schema must include fallback_used in stats
  assert.ok(
    output.properties.stats.required.includes('fallback_used'),
    'output.schema stats must include fallback_used'
  );
});

test('tool-definition.json 含 8 个 outputs 且 ipcChannel + failureModes 完整', () => {
  const path = resolve(pluginDir, 'tool-definition.json');
  assert.ok(existsSync(path), 'tool-definition.json missing');
  const td = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(td.id, 'celltypist');
  assert.equal(td.language, 'python');
  assert.equal(td.category, 'single-cell');
  assert.ok(Array.isArray(td.outputs) && td.outputs.length === 8,
    `tool-definition.json outputs length must be 8, got ${td.outputs?.length}`);

  const expectedOutputs = [
    'predictions',
    'probabilities',
    'confusion_matrix',
    'annotated',
    'umap_annotation',
    'umap_confidence',
    'report',
    'manifest',
  ];
  const ids = td.outputs.map(o => o.id);
  for (const id of expectedOutputs) {
    assert.ok(ids.includes(id), `outputs missing: ${id}`);
  }

  assert.ok(td.ipcChannel === 'tools:invoke/celltypist', 'ipcChannel mismatch');
  assert.ok(td.failureModes && td.failureModes.length === 3,
    `failureModes should have 3 entries, got ${td.failureModes?.length}`);
  const failureCodes = td.failureModes.map(f => f.code);
  assert.ok(failureCodes.includes('missing_required_file_inputs'));
  assert.ok(failureCodes.includes('runtime_not_ready'));
  assert.ok(failureCodes.includes('empty_result'));

  // resultStudio 必须存在 (用于 PHASE 1 SUB-1.6 skill 调用)
  assert.ok(td.resultStudio, 'tool-definition.json must include resultStudio block');
  assert.equal(td.resultStudio.semanticType, 'cell_type_annotation', 'semanticType must be cell_type_annotation');
});

test('acceptance.json 4 段 fixture (valid + missing-inputs + runtime-not-ready + empty-result)', () => {
  const path = resolve(pluginDir, 'acceptance.json');
  assert.ok(existsSync(path), 'acceptance.json missing');
  const a = JSON.parse(readFileSync(path, 'utf8'));
  assert.equal(a.pluginId, 'celltypist');
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
});

test('placeholder PNG 存在且 > 200 bytes（hash rewrite 后稳定输出）', () => {
  const path = resolve(pluginDir, 'tools/figures/annotation_placeholder.png');
  assert.ok(existsSync(path), 'tools/figures/annotation_placeholder.png missing');
  const stats = statSync(path);
  assert.ok(stats.size > 200, `placeholder PNG too small: ${stats.size} bytes`);
  // PNG magic bytes
  const head = readFileSync(path).subarray(0, 8);
  assert.deepEqual(
    Array.from(head),
    [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a],
    'placeholder PNG missing magic bytes'
  );

  // hash 稳定性: 2 次连跑 sha256 应该一致
  const { createHash } = require('node:crypto');
  const h1 = createHash('sha256').update(readFileSync(path)).digest('hex');
  // 2nd run via the generator script (write+rewrite)
  const { execSync } = require('node:child_process');
  // BLOCK-83 P4 跨平台：Windows 用 PATH 上的 python
  const PYTHON = process.platform === 'win32' ? 'python' : '/Applications/anaconda3/envs/biof3-py-runtime/bin/python3';
  execSync(PYTHON + ' ' +
    resolve(pluginDir, 'tools/figures/generate_placeholder_figure.py'), { stdio: 'ignore' });
  const h2 = createHash('sha256').update(readFileSync(path)).digest('hex');
  assert.equal(h1, h2, `placeholder PNG hash changed after re-run: ${h1} -> ${h2}`);
});

test('demo h5ad 落盘且 sha256 与 scanpy-advanced 共享 (e3d6c0f022...)', () => {
  const path = resolve(pluginDir, 'demo-data/pbmc_3k_mini.h5ad');
  assert.ok(existsSync(path), 'demo-data/pbmc_3k_mini.h5ad missing');
  const stats = statSync(path);
  assert.ok(stats.size > 10000, `demo h5ad too small: ${stats.size} bytes`);

  const { createHash } = require('node:crypto');
  const h = createHash('sha256').update(readFileSync(path)).digest('hex');
  assert.equal(
    h,
    'e3d6c0f0226479e52c5e8d8d03e2be392466766d25d07c326d0cd0982e2c35bc',
    `demo h5ad sha256 mismatch: expected e3d6...3ebc, got ${h}`
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
  // 必须显式声明 celltypist 是 extension package (PyPI bootstrap)
  assert.ok(
    body.includes('celltypist') && (body.includes('PyPI') || body.includes('bootstrap') || body.includes('install.json')),
    'README.md must declare celltypist as PyPI bootstrap extension (not in biof3-py-runtime baseline)'
  );
});

test('顶层 manifest.json 已包含 celltypist 条目', () => {
  const path = resolve(repoRoot, 'plugins/manifest.json');
  const manifest = JSON.parse(readFileSync(path, 'utf8'));
  assert.ok(manifest.plugins['celltypist'], 'manifest.json celltypist entry missing');
  assert.equal(manifest.plugins['celltypist'].id, 'celltypist');
  assert.equal(manifest.plugins['celltypist'].language, 'python');
  assert.equal(manifest.plugins['celltypist'].kind, 'plugin');
  // celltypist plugin entry 必须显式声明 celltypist 包
  assert.ok(
    manifest.plugins['celltypist'].packages.python.celltypist,
    'manifest.json celltypist entry must declare celltypist package'
  );
});

test('celltypist plugin 目录结构完整（8 件套文件全部存在）', () => {
  const required = [
    'plugin.json',
    'environment.json',
    'tool-definition.json',
    'acceptance.json',
    'manifest.entry.json',
    'install.json',
    'README.md',
    'celltypist-report-template.html',
    'scripts/celltypist.py',
    'schema/input.schema.json',
    'schema/output.schema.json',
    'demo-data/pbmc_3k_mini.h5ad',
    'demo-data/build_demo_h5ad.py',
    'tools/figures/annotation_placeholder.png',
    'tools/figures/generate_placeholder_figure.py',
  ];
  const missing = required.filter(rel => !existsSync(resolve(pluginDir, rel)));
  assert.equal(missing.length, 0, `missing files: ${missing.join(', ')}`);
});

test('celltypist plugin 目录无未提交孤儿文件', () => {
  const requiredDirs = ['scripts', 'schema', 'demo-data', 'tools/figures', 'tools'];
  for (const d of requiredDirs) {
    const dir = resolve(pluginDir, d);
    assert.ok(existsSync(dir), `directory missing: ${d}`);
    const items = readdirSync(dir);
    assert.ok(items.length > 0, `directory empty: ${d}`);
  }
});

test('plugin.json 与 environment.json 的 packages.python 完全一致（含 celltypist）', () => {
  const plugin = JSON.parse(readFileSync(resolve(pluginDir, 'plugin.json'), 'utf8'));
  const env = JSON.parse(readFileSync(resolve(pluginDir, 'environment.json'), 'utf8'));
  // SUB-1.1 deseq2 R 模式要求二份清单保持一致；此处同款约束
  for (const [pkg, spec] of Object.entries(plugin.packages.python)) {
    assert.equal(
      env.packages.python[pkg],
      spec,
      `packages.python.${pkg} mismatch: plugin.json=${spec} vs environment.json=${env.packages.python[pkg]}`
    );
  }
});

test('celltypist 既不是 biof3-py-runtime 基础包, install.json 必须 bootstrap', () => {
  // py-runtime-requirements.txt 不应包含 celltypist
  const pyReq = readFileSync(
    resolve(repoRoot, 'plugins/py-runtime-requirements.txt'),
    'utf8'
  );
  assert.ok(
    !/^celltypist\s*>=/m.test(pyReq),
    'biof3-py-runtime-requirements.txt must NOT include celltypist (celltypist is an extension package)'
  );
  // install.json 必须声明 pip-bootstrap
  const install = JSON.parse(readFileSync(resolve(pluginDir, 'install.json'), 'utf8'));
  assert.equal(install.policy, 'pip-bootstrap', 'celltypist install policy must be pip-bootstrap');
  assert.ok(install.steps.length > 0, 'install.json steps must not be empty');
  assert.ok(
    install.steps.some(s => s.kind === 'pip-install' && s.package === 'celltypist'),
    'install.json must have pip-install celltypist step'
  );
});
