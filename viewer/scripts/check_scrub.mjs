// Scrub the time slider through every saved step of the loaded run; fail on any page error
// or if the status line does not follow the slider.
import { chromium } from 'playwright';
const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1200, height: 800 } });
const errors = [];
page.on('pageerror', (e) => errors.push(String(e)));
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
await page.goto(process.argv[2] ?? 'http://localhost:5173/');
await page.waitForTimeout(1500);
const slider = page.locator('input.pd-scrub');
const max = Number(await slider.getAttribute('max'));
let bad = 0;
for (let i = 0; i <= max; i++) {
  await slider.fill(String(i));
  await page.waitForTimeout(60);
  const label = await page.locator('span', { hasText: /^t = / }).first().textContent();
  const t = Number(label.replace('t = ', ''));
  await page.waitForFunction((tt) => document.body.innerText.includes(`t = ${tt} Myr`), t, { timeout: 3000 })
    .catch(() => { bad++; });
}
await browser.close();
console.log(`${max + 1} steps scrubbed, ${bad} status mismatches, ${errors.length} errors`);
if (errors.length || bad) { console.error(errors.join('\n')); process.exit(1); }
