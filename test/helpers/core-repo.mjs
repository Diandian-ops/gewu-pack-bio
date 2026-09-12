import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const PACK_ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..');

/** 定位核心仓（兄弟目录 biof3-desktop；可用 BIOF3_CORE_REPO 覆盖）。 */
export function findCoreRepo() {
  if (process.env.BIOF3_CORE_REPO && fs.existsSync(process.env.BIOF3_CORE_REPO)) {
    return process.env.BIOF3_CORE_REPO;
  }
  const sibling = path.resolve(PACK_ROOT, '..', 'biof3-desktop');
  if (fs.existsSync(path.join(sibling, 'package.json'))) return sibling;
  return null;
}

export const PACK_ROOT_EXPORT = PACK_ROOT;
