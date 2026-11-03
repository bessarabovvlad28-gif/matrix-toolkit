const EPS = 1e-10;
const $ = (id) => document.getElementById(id);
const output = $('result');

function parse(text) {
  const rows = text.trim().split(/[\n;]+/).filter(Boolean).map(line => line.trim().split(/[\s,]+/).map(Number));
  if (!rows.length || rows.some(row => !row.length || row.some(value => !Number.isFinite(value)))) throw Error('Введите матрицу из чисел. Разделяйте строки переносом, числа — пробелом.');
  if (rows.some(row => row.length !== rows[0].length)) throw Error('Во всех строках матрицы должно быть одинаковое количество чисел.');
  return rows;
}
const fmt = n => Math.abs(n) < EPS ? '0' : String(Number(n.toFixed(6)));
const showMatrix = m => m.map(row => row.map(fmt).join('   ')).join('\n');
const same = (a,b) => a.length === b.length && a[0].length === b[0].length;
function transpose(a) { return a[0].map((_,j) => a.map(row => row[j])); }
function add(a,b,sign=1) { if (!same(a,b)) throw Error('Для этой операции размеры A и B должны совпадать.'); return a.map((r,i) => r.map((v,j) => v + sign*b[i][j])); }
function multiply(a,b) { if (a[0].length !== b.length) throw Error('Число столбцов A должно совпадать с числом строк B.'); return a.map(row => b[0].map((_,j) => row.reduce((sum,v,k) => sum + v*b[k][j],0))); }
function rref(input) {
  const a = input.map(row => [...row]); let lead = 0;
  for (let r=0; r<a.length && lead<a[0].length; r++) {
    let p=r; while (p<a.length && Math.abs(a[p][lead])<EPS) p++;
    if (p===a.length) { lead++; r--; continue; }
    [a[r],a[p]]=[a[p],a[r]]; const pivot=a[r][lead]; a[r]=a[r].map(v=>v/pivot);
    for (let i=0;i<a.length;i++) if(i!==r) { const f=a[i][lead]; a[i]=a[i].map((v,j)=>v-f*a[r][j]); }
    lead++;
  }
  return a.map(row=>row.map(v=>Math.abs(v)<EPS?0:v));
}
function det(a) { if(a.length!==a[0].length) throw Error('Определитель существует только у квадратной матрицы.'); const m=a.map(r=>[...r]); let d=1;
  for(let i=0;i<m.length;i++){ let p=i; for(let r=i+1;r<m.length;r++) if(Math.abs(m[r][i])>Math.abs(m[p][i]))p=r; if(Math.abs(m[p][i])<EPS)return 0; if(p!==i){[m[i],m[p]]=[m[p],m[i]];d=-d;} const x=m[i][i];d*=x;for(let r=i+1;r<m.length;r++){const f=m[r][i]/x;for(let c=i+1;c<m.length;c++)m[r][c]-=f*m[i][c];}}
  return d;
}
function inverse(a) { if(a.length!==a[0].length) throw Error('Обратная матрица существует только для квадратной матрицы.'); const n=a.length,m=a.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>+(i===j))]);
  for(let i=0;i<n;i++){let p=i;for(let r=i+1;r<n;r++)if(Math.abs(m[r][i])>Math.abs(m[p][i]))p=r;if(Math.abs(m[p][i])<EPS)throw Error('Матрица вырожденная: обратной матрицы нет.');[m[i],m[p]]=[m[p],m[i]];const x=m[i][i];m[i]=m[i].map(v=>v/x);for(let r=0;r<n;r++)if(r!==i){const f=m[r][i];m[r]=m[r].map((v,j)=>v-f*m[i][j]);}}
  return m.map(r=>r.slice(n));
}
function solve(a,b) { if(a.length!==a[0].length) throw Error('Для Ax = b матрица A должна быть квадратной.'); if(b[0].length!==1||b.length!==a.length)throw Error('Для Ax = b введите в B один столбец, по одному числу на строку (столько же строк, сколько у A).'); const rr=rref(a.map((r,i)=>[...r,b[i][0]])); const n=a.length; if(rr.some(r=>r.slice(0,n).every(v=>Math.abs(v)<EPS)&&Math.abs(r[n])>EPS))return 'Система не имеет решений.'; if(rr.filter(r=>r.slice(0,n).some(v=>Math.abs(v)>EPS)).length<n)return 'Система имеет бесконечно много решений.'; return 'x = ('+rr.map(r=>fmt(r[n])).join(', ')+')'; }
function run(event) {
  try {
    const a=parse($('matrixA').value), b=parse($('matrixB').value), op=event.currentTarget.dataset.op; let value;
    switch(op){
      case 'add': value=showMatrix(add(a,b)); break; case 'subtract': value=showMatrix(add(a,b,-1)); break; case 'multiply': value=showMatrix(multiply(a,b)); break;
      case 'transposeA': value=showMatrix(transpose(a)); break; case 'transposeB': value=showMatrix(transpose(b)); break;
      case 'detA': value=`det(A) = ${fmt(det(a))}`; break; case 'detB': value=`det(B) = ${fmt(det(b))}`; break;
      case 'inverseA': value=showMatrix(inverse(a)); break; case 'inverseB': value=showMatrix(inverse(b)); break;
      case 'rankA': value=`rank(A) = ${rref(a).filter(r=>r.some(v=>Math.abs(v)>EPS)).length}`; break;
      case 'traceA': if(a.length!==a[0].length)throw Error('След вычисляется только у квадратной матрицы.'); value=`trace(A) = ${fmt(a.reduce((s,r,i)=>s+r[i],0))}`; break;
      case 'rrefA': value=showMatrix(rref(a)); break; case 'solve': value=solve(a,b); break;
      default: value='Неизвестная операция.';
    }
    output.textContent=value;
  } catch(error) { output.textContent=error.message; }
}
document.querySelectorAll('[data-op]').forEach(button=>button.addEventListener('click',run));
$('sample').addEventListener('click',()=>{$('matrixA').value='2 1 0\n1 3 1\n0 1 2';$('matrixB').value='1 0 0\n0 1 0\n0 0 1';output.textContent='Пример загружен. Выберите операцию.';});
