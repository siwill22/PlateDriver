// Switch projection with the real <select> in the panel (Globe -> Robinson -> Plate Carree
// -> Globe), then scrub; screenshot each; fail on page errors.
import { chromium } from 'playwright';
const out = process.argv[2] ?? 'shots/switch';
const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1200, height: 800 } });
const errors = [];
page.on('pageerror', (e) => errors.push(String(e)));
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
await page.goto('http://localhost:5173/?t=250');
await page.waitForTimeout(1500);
const select = page.locator('.lil-gui select').first();
for (const label of ['Robinson', 'Plate Carrée', 'Globe', 'Robinson']) {
  await select.selectOption({ label });
  await page.waitForTimeout(700);
  await page.screenshot({ path: `${out}_${label.replace(/\W/g, '')}.png` });
}
await page.locator('input.pd-scrub').fill('50');
await page.waitForTimeout(700);
await page.screenshot({ path: `${out}_Robinson_t500.png` });
await browser.close();
console.log(`${errors.length} errors`);
if (errors.length) { console.error(errors.join('\n')); process.exit(1); }
