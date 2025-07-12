// -- offline template resolver (user-provided shim) -----------------------------
const LOCAL_TPL = Object.fromEntries(
  [...document.querySelectorAll('script[type="text/template"]')]
    .map(t => [t.id, t.innerHTML])
);

async function fetchTpl(path) {
  const m = /src\/templates\/([^/]+\.tpl)$/.exec(path);
  if (m && LOCAL_TPL[m[1]]) {
    return new Response(LOCAL_TPL[m[1]], { status: 200 });
  }
  // Fallback to actual fetch for any non-template request
  return fetch(path);
}


// https://stackoverflow.com/a/73891404
async function replaceAsync(string, regexp, replacerFunction) {
    const replacements = await Promise.all(
        Array.from(string.matchAll(regexp), replacerFunction));
    let i = 0;

    return string.replace(regexp, () => replacements[i++]);
}


function compileTemplate(template) {
    return replaceAsync(template, /\{\{(\S+)\}\}/g, async match => {
        // Use the fetchTpl shim instead of the original fetch
        return await fetchTpl(`src/templates/${match[1]}`)
            .then(response => response.text())
            .then(template => {
                return compileTemplate(template);
            });
        }
    );
}

const main = document.querySelector('main');

compileTemplate(main.innerHTML).then((result) => {
    main.innerHTML = result;
    import('./main.js');
});

