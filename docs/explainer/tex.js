// LaTeX to MathML with Temml: reads {name: [tex, display]} on stdin, writes {name: mathml}.
const temml = require('temml');
const t = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const out = {};
for (const [k, [tex, display]] of Object.entries(t))
  out[k] = temml.renderToString(tex, {displayMode: display, throwOnError: true});
process.stdout.write(JSON.stringify(out));
