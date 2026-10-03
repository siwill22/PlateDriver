// Switch runs with the panel's Run dropdown; fail on any page error or if the status line
// does not show the chosen run's label.
import { chromium } from 'playwright';
const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1200, height: 800 } });
const errors = [];
page.on('pageerror', (e) => errors.push(String(e)));
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
await page.goto(process.argv[2] ?? 'http://localhost:5173/');
await page.waitForTimeout(1500);
const select = page.locator('.lil-gui select').last();
const labels = await select.locator('option').allTextContents();
let bad = 0;
for (const label of labels) {
  await select.selectOption({ label });
  await page.waitForFunction((l) => document.body.innerText.includes(l), label, { timeout: 5000 })
    .catch(() => { bad++; console.error(`status never showed ${label}`); });
}
await browser.close();
console.log(`${labels.length} runs switched, ${bad} mismatches, ${errors.length} errors`);
if (errors.length || bad) { console.error(errors.join('\n')); process.exit(1); }
