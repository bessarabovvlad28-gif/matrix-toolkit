const EPS = 1e-10;
const $ = (id) => document.getElementById(id);
const output = $('result');
const stepsOutput = $('steps');

function parse(text) {
  const rows = text.trim().split(/[\n;]+/).filter(Boolean).map(line => line.trim().split(/[\s,]+/).map(Number));
  if (!rows.length || rows.some(row => !row.length || row.some(value => !Number.isFinite(value)))) {
    throw Error('Введите матрицу из чисел. Разделяйте строки переносом, числа — пробелом.');
  }
  if (rows.some(row => row.length !== rows[0].length)) {
    throw Error('Во всех строках матрицы должно быть одинаковое количество чисел.');
  }
  return rows;
}

const fmt = n => Math.abs(n) < EPS ? '0' : String(Number(n.toFixed(6)));
const showMatrix = m => m.map(row => row.map(fmt).join('   ')).join('\n');
const same = (a,b) => a.length === b.length && a[0].length === b[0].length;
function transpose(a) { return a[0].map((_,j) => a.map(row => row[j])); }
function add(a,b,sign=1) {
  if (!same(a,b)) throw Error('Для этой операции размеры A и B должны совпадать.');
  return a.map((r,i) => r.map((v,j) => v + sign*b[i][j]));
}
function multiply(a,b) {
  if (a[0].length !== b.length) throw Error('Число столбцов A должно совпадать с числом строк B.');
  return a.map(row => b[0].map((_,j) => row.reduce((sum,v,k) => sum + v*b[k][j],0)));
}
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
function det(a) {
  if(a.length!==a[0].length) throw Error('Определитель существует только у квадратной матрицы.');
  const m=a.map(r=>[...r]); let result=1; let sign=1;
  for(let i=0;i<m.length;i++) {
    let p=i; for(let r=i+1;r<m.length;r++) if(Math.abs(m[r][i])>Math.abs(m[p][i])) p=r;
    if(Math.abs(m[p][i])<EPS) return 0;
    if(p!==i) { [m[i],m[p]]=[m[p],m[i]]; sign*=-1; }
    const pivot=m[i][i]; result*=pivot;
    for(let r=i+1;r<m.length;r++) { const f=m[r][i]/pivot; for(let c=i+1;c<m.length;c++) m[r][c]-=f*m[i][c]; }
  }
  return result*sign;
}
function inverse(a) {
  if(a.length!==a[0].length) throw Error('Обратная матрица существует только для квадратной матрицы.');
  const n=a.length,m=a.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>+(i===j))]);
  for(let i=0;i<n;i++) {
    let p=i; for(let r=i+1;r<n;r++) if(Math.abs(m[r][i])>Math.abs(m[p][i])) p=r;
    if(Math.abs(m[p][i])<EPS) throw Error('Матрица вырожденная: обратной матрицы нет.');
    [m[i],m[p]]=[m[p],m[i]]; const x=m[i][i]; m[i]=m[i].map(v=>v/x);
    for(let r=0;r<n;r++) if(r!==i) { const f=m[r][i]; m[r]=m[r].map((v,j)=>v-f*m[i][j]); }
  }
  return m.map(r=>r.slice(n));
}
function trace(a) {
  if(a.length!==a[0].length) throw Error('След вычисляется только у квадратной матрицы.');
  return a.reduce((s,r,i)=>s+r[i],0);
}

function parseSystem() {
  const A=parse($('matrixA').value), rhs=parse($('matrixB').value);
  if(A.length!==A[0].length) throw Error('Для методов Крамера и Гаусса матрица A должна быть квадратной.');
  let b;
  if(rhs[0].length===1) b=rhs.map(row=>row[0]);
  else if(rhs.length===1 && rhs[0].length===A.length) b=[...rhs[0]];
  else throw Error('Введите b одним столбцом (по одному числу на строку) или одной строкой.');
  if(b.length!==A.length) throw Error(`В векторе b должно быть ${A.length} чисел — по одному на каждую строку A.`);
  return {A,b};
}

function cramer(A,b) {
  const D=det(A);
  const log=[`Метод Крамера для системы Ax = b`,`\n1. Матрица коэффициентов A:\n${showMatrix(A)}`,`\n2. Определитель матрицы A: D = det(A) = ${fmt(D)}.`];
  if(Math.abs(D)<EPS) {
    log.push('Так как D = 0, метод Крамера неприменим: система не имеет единственного решения.');
    const rankA=rref(A).filter(row=>row.some(v=>Math.abs(v)>EPS)).length;
    const augmented=A.map((row,i)=>[...row,b[i]]);
    const rankAug=rref(augmented).filter(row=>row.some(v=>Math.abs(v)>EPS)).length;
    const answer=rankA<rankAug ? 'Система несовместна (решений нет).' : 'Система совместна, но единственного решения нет (решений бесконечно много).';
    return {answer,steps:log.join('\n')};
  }
  const x=[];
  for(let j=0;j<A.length;j++) {
    const Aj=A.map((row,i)=>row.map((v,c)=>c===j?b[i]:v));
    const Dj=det(Aj); x.push(Dj/D);
    log.push(`\n3.${j+1}. Заменяем ${j+1}-й столбец A на b:\nA${j+1} =\n${showMatrix(Aj)}\nD${j+1} = det(A${j+1}) = ${fmt(Dj)}\nx${j+1} = D${j+1} / D = ${fmt(Dj)} / ${fmt(D)} = ${fmt(x[j])}.`);
  }
  log.push(`\nОтвет: x = (${x.map(fmt).join(', ')}).`);
  return {answer:`Решение: x = (${x.map(fmt).join(', ')})`,steps:log.join('\n')};
}

function gauss(A,b) {
  if(A.length!==A[0].length) throw Error('Метод Гаусса сейчас рассчитан на квадратную матрицу A.');
  const n=A.length, m=A.map((row,i)=>[...row,b[i]]), log=['Метод Гаусса: составляем расширенную матрицу [A | b].',showMatrix(m)];
  let pivotRow=0;
  for(let col=0;col<n && pivotRow<n;col++) {
    let pivot=pivotRow;
    for(let r=pivotRow+1;r<n;r++) if(Math.abs(m[r][col])>Math.abs(m[pivot][col])) pivot=r;
    if(Math.abs(m[pivot][col])<EPS) { log.push(`\nВ столбце ${col+1} нет ведущего элемента; переходим к следующему столбцу.`); continue; }
    if(pivot!==pivotRow) { [m[pivotRow],m[pivot]]=[m[pivot],m[pivotRow]]; log.push(`\nМеняем местами строки R${pivotRow+1} и R${pivot+1}:\n${showMatrix(m)}`); }
    const pivotValue=m[pivotRow][col];
    log.push(`\nВедущий элемент в строке ${pivotRow+1}, столбце ${col+1}: ${fmt(pivotValue)}.`);
    for(let r=pivotRow+1;r<n;r++) {
      const factor=m[r][col]/pivotValue;
      if(Math.abs(factor)<EPS) continue;
      for(let c=col;c<=n;c++) m[r][c]-=factor*m[pivotRow][c];
      m[r]=m[r].map(v=>Math.abs(v)<EPS?0:v);
      log.push(`R${r+1} ← R${r+1} − (${fmt(factor)})·R${pivotRow+1}\n${showMatrix(m)}`);
    }
    pivotRow++;
  }
  for(let r=0;r<n;r++) {
    if(m[r].slice(0,n).every(v=>Math.abs(v)<EPS) && Math.abs(m[r][n])>=EPS) {
      log.push('\nПолучено равенство 0 = ненулевое число, поэтому система несовместна.');
      return {answer:'Система не имеет решений.',steps:log.join('\n')};
    }
  }
  const pivots=[];
  for(let r=0;r<n;r++) { const c=m[r].slice(0,n).findIndex(v=>Math.abs(v)>EPS); if(c>=0) pivots.push([r,c]); }
  if(pivots.length<n) {
    log.push('\nВедущих переменных меньше, чем неизвестных. Система имеет бесконечно много решений.');
    return {answer:'Система имеет бесконечно много решений.',steps:log.join('\n')};
  }
  const x=Array(n).fill(0);
  log.push('\nОбратный ход: выражаем неизвестные, начиная с последней строки.');
  for(let r=n-1;r>=0;r--) {
    let value=m[r][n];
    for(let c=r+1;c<n;c++) value-=m[r][c]*x[c];
    x[r]=value/m[r][r];
    log.push(`x${r+1} = ${fmt(x[r])}`);
  }
  log.push(`\nОтвет: x = (${x.map(fmt).join(', ')}).`);
  return {answer:`Решение: x = (${x.map(fmt).join(', ')})`,steps:log.join('\n')};
}

function run(event) {
  try {
    const a=parse($('matrixA').value), b=parse($('matrixB').value), op=event.currentTarget.dataset.op;
    let value='';
    if(op==='cramer'||op==='gauss') {
      const system=parseSystem();
      const result=op==='cramer'?cramer(system.A,system.b):gauss(system.A,system.b);
      output.textContent=result.answer; stepsOutput.textContent=result.steps; return;
    }
    stepsOutput.textContent='Подробный ход решения пока показывается для методов Крамера и Гаусса.';
    switch(op) {
      case 'add': value=showMatrix(add(a,b)); break;
      case 'subtract': value=showMatrix(add(a,b,-1)); break;
      case 'multiply': value=showMatrix(multiply(a,b)); break;
      case 'transposeA': value=showMatrix(transpose(a)); break;
      case 'transposeB': value=showMatrix(transpose(b)); break;
      case 'detA': value=`det(A) = ${fmt(det(a))}`; break;
      case 'detB': value=`det(B) = ${fmt(det(b))}`; break;
      case 'inverseA': value=showMatrix(inverse(a)); break;
      case 'inverseB': value=showMatrix(inverse(b)); break;
      case 'rankA': value=`rank(A) = ${rref(a).filter(r=>r.some(v=>Math.abs(v)>EPS)).length}`; break;
      case 'traceA': value=`trace(A) = ${fmt(trace(a))}`; break;
      case 'rrefA': value=showMatrix(rref(a)); break;
      default: value='Неизвестная операция.';
    }
    output.textContent=value;
  } catch(error) { output.textContent=error.message; stepsOutput.textContent=''; }
}

document.querySelectorAll('[data-op]').forEach(button=>button.addEventListener('click',run));
$('sample').addEventListener('click',()=>{
  $('matrixA').value='2 1 0\n1 3 1\n0 1 2';
  $('matrixB').value='1\n2\n3';
  output.textContent='Пример загружен. Выберите метод Крамера или Гаусса.';
  stepsOutput.textContent='';
});
