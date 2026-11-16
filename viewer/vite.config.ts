import { defineConfig } from 'vite';

// Local only for now (see docs/viewer_plan.md): `npm run dev` on localhost:5173.
export default defineConfig({
  base: process.env.VITE_BASE ?? '/',
  server: { port: 5173 },
  build: { target: 'es2022', chunkSizeWarningLimit: 700 },
});
