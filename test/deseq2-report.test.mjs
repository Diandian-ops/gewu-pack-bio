// test/deseq2-report.test.mjs — #49 回归：报告模板填充 + 平台安全模板解析
//
// 锁定两点：
//  1. deseq2-report-template.html 必须包含 7 张图的 <img> 锚点 URL（deseq2.R 的 gsub
//     按这些 URL 替换为 base64）以及 9 个摘要占位符 {{...}}，否则报告仍是空壳。
//  2. deseq2.R 不得再硬编码 POSIX 路径 /opt/biof3-biof3-server/...（Windows/macOS app
//     bundle 均不存在）；改用宿主注入的 BIOF3_BUILT_IN_PLUGINS_DIR 兜底，并实际做
//     摘要 token 替换。

import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..'); // BLOCK-83 P4：包根

const templatePath = require('path').join(
  repoRoot, 'plugins/deseq2/deseq2-report-template.html'
);
const scriptPath = require('path').join(
  repoRoot, 'plugins/deseq2/scripts/deseq2.R'
);

const template = readFileSync(templatePath, 'utf8');
const script = readFileSync(scriptPath, 'utf8');

// ── 1. 模板锚点 + 占位符 ──────────────────────────────────────
const EXPECTED_IMG_URLS = [
  'library_size.png',
  'dispersion.png',
  'pca.png',
  'sample_distance.png',
  'volcano.png',
  'ma_plot.png',
  'heatmap_top50.png', // 注意：R 脚本 img_files 用此 key（映射到 heatmap.png 文件）
].map((n) => `https://biof3.com/api/r/tools/demo/deseq2/${n}`);

test('report template embeds all 7 plot <img> anchors', () => {
  for (const url of EXPECTED_IMG_URLS) {
    assert.ok(
      template.includes(`src="${url}"`),
      `template missing <img> anchor for ${url}`
    );
  }
});

const SUMMARY_TOKENS = [
  '{{TOTAL_GENES}}', '{{SIG_GENES}}', '{{UP}}', '{{DOWN}}',
  '{{DESIGN}}', '{{REF}}', '{{TREAT}}', '{{PADJ}}', '{{LFC}}',
];
test('report template declares all 9 summary tokens', () => {
  for (const t of SUMMARY_TOKENS) {
    assert.ok(template.includes(t), `template missing summary token ${t}`);
  }
});

test('report template declares the DEG results-table token', () => {
  assert.ok(
    template.includes('{{DEG_TABLE}}'),
    'template missing {{DEG_TABLE}} token for the gene-level results table'
  );
});

// ── 2. R 脚本：无 POSIX 硬编码 + 平台安全兜底 + token 替换 ──
test('deseq2.R does NOT hardcode POSIX /opt path', () => {
  assert.ok(
    !script.includes('/opt/biof3-biof3-server'),
    'deseq2.R must not hardcode the POSIX /opt/biof3-biof3-server path (broken on Windows/macOS bundle)'
  );
});

test('deseq2.R uses platform-safe BIOF3_BUILT_IN_PLUGINS_DIR fallback', () => {
  assert.ok(
    script.includes('BIOF3_BUILT_IN_PLUGINS_DIR'),
    'deseq2.R should resolve missing template via BIOF3_BUILT_IN_PLUGINS_DIR env (set by main.js)'
  );
  assert.ok(
    script.includes('file.path(built_in_dir, "deseq2", "deseq2-report-template.html")'),
    'fallback must use cross-platform file.path, not a hardcoded POSIX string'
  );
});

test('deseq2.R performs summary token substitution', () => {
  for (const t of SUMMARY_TOKENS) {
    assert.ok(
      script.includes(`gsub("${t}"`),
      `deseq2.R should replace summary token ${t}`
    );
  }
});

test('deseq2.R injects the DEG results table', () => {
  assert.ok(
    script.includes('gsub("{{DEG_TABLE}}"'),
    'deseq2.R should replace the {{DEG_TABLE}} token with <tr> rows built from sig'
  );
});
