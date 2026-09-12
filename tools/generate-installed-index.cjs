#!/usr/bin/env node
'use strict';

/**
 * 生成包的插件安装索引 `plugins/installed.json`。
 *
 * 核心的 installed 源布局要求：`<root>/installed.json`（数组）+ `<root>/<id>/plugin.json`。
 * 本脚本从 plugins/ 下的各插件目录汇总索引，字段与核心一致：
 *   { id, version, installedAt, platform, language, type }
 *
 * 用法：node tools/generate-installed-index.cjs
 */

const fs = require('fs');
const path = require('path');

const PACK_ROOT = path.resolve(__dirname, '..');
const PLUGINS_DIR = path.join(PACK_ROOT, 'plugins');
const OUT = path.join(PLUGINS_DIR, 'installed.json');

if (!fs.existsSync(PLUGINS_DIR)) {
  console.error('plugins/ 目录不存在: ' + PLUGINS_DIR);
  process.exit(2);
}

const entries = [];
for (const entry of fs.readdirSync(PLUGINS_DIR, { withFileTypes: true })) {
  if (!entry.isDirectory()) continue;
  const pluginJson = path.join(PLUGINS_DIR, entry.name, 'plugin.json');
  if (!fs.existsSync(pluginJson)) continue;
  const plugin = JSON.parse(fs.readFileSync(pluginJson, 'utf8'));
  entries.push({
    id: plugin.id || entry.name,
    version: plugin.version || '0.0.0',
    installedAt: new Date().toISOString(),
    platform: process.platform === 'win32' ? 'windows-x64' : process.platform === 'darwin' ? 'mac-arm64' : 'linux-x64',
    language: plugin.language || null,
    type: plugin.type || 'analysis-plugin',
  });
}

fs.writeFileSync(OUT, JSON.stringify(entries, null, 2) + '\n', 'utf8');
console.log('[pack] 写入 ' + entries.length + ' 条插件索引 → plugins/installed.json');
for (const e of entries) console.log('  ' + e.id + ' @ ' + e.version + ' (' + e.type + ')');
