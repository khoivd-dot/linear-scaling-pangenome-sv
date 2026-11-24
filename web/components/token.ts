// Resolve a CSS custom property to rgb(): Cytoscape cannot parse oklch, and canvas
// drawing should use the same palette as globals.css.
export function token(name: string): string {
  const value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const ctx = document.createElement("canvas").getContext("2d", { willReadFrequently: true });
  if (!ctx || !value) return value;
  ctx.fillStyle = value;
  ctx.fillRect(0, 0, 1, 1);
  const d = ctx.getImageData(0, 0, 1, 1).data;
  return `rgb(${d[0]}, ${d[1]}, ${d[2]})`;
}
