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
const showMatrix = matrix => matrix.map(row => row.map(fmt).join('   ')).join('\n');
const sameSize = (a,b) => a.length === b.length && a[0].length === b[0].length;
const clone = a => a.map(row => [...row]);

function parseSystem() {
  const A = parse($('matrixA').value);
  const rhs = parse($('matrixB').value);
  if (A.length !== A[0].length) throw Error('Для методов Крамера и Гаусса матрица A должна быть квадратной.');
  let b;
  if (rhs[0].length === 1) b = rhs.map(row => row[0]);
  else if (rhs.length === 1 && rhs[0].length === A.length) b = [...rhs[0]];
  else throw Error('Введите b одним столбцом или одной строкой.');
  if (b.length !== A.length) throw Error(`В векторе b должно быть ${A.length} чисел.`);
  return { A, b };
}

function transpose(a) {
  return a[0].map((_, j) => a.map(row => row[j]));
}

function add(a, b, sign = 1) {
  if (!sameSize(a, b)) throw Error('Размеры матриц должны совпадать.');
  return a.map((row, i) => row.map((v, j) => v + sign * b[i][j]));
}

function multiply(a, b) {
  if (a[0].length !== b.length) throw Error('Число столбцов A должно совпадать с числом строк B.');
  return a.map(row => b[0].map((_, j) => row.reduce((sum, v, k) => sum + v * b[k][j], 0)));
}

function determinantSteps(input) {
  if (input.length !== input[0].length) throw Error('Определитель рассчитывается только для квадратной матрицы.');
  const a = clone(input), n = a.length, log = [`Исходная матрица:\n${showMatrix(a)}`];
  let sign = 1, product = 1;
  for (let i = 0; i < n; i++) {
    let pivotRow = i;
    for (let r = i + 1; r < n; r++) if (Math.abs(a[r][i]) > Math.abs(a[pivotRow][i])) pivotRow = r;
    if (Math.abs(a[pivotRow][i]) < EPS) {
      log.push(`В столбце ${i + 1} ведущий элемент равен нулю: определитель равен 0.`);
      return { value: 0, steps: log.join('\n\n') };
    }
    if (pivotRow !== i) {
      [a[i], a[pivotRow]] = [a[pivotRow], a[i]];
      sign *= -1;
      log.push(`Меняем строки ${i + 1} и ${pivotRow + 1}; знак определителя меняется.\n${showMatrix(a)}`);
    }
    const pivot = a[i][i];
    product *= pivot;
    log.push(`Ведущий элемент a${i + 1}${i + 1} = ${fmt(pivot)}. Обнуляем элементы ниже него.`);
    for (let r = i + 1; r < n; r++) {
      const factor = a[r][i] / pivot;
      if (Math.abs(factor) < EPS) continue;
      for (let c = i; c < n; c++) a[r][c] -= factor * a[i][c];
      a[r] = a[r].map(v => Math.abs(v) < EPS ? 0 : v);
      log.push(`R${r + 1} ← R${r + 1} − (${fmt(factor)})·R${i + 1}:\n${showMatrix(a)}`);
    }
  }
  const value = sign * product;
  log.push(`Определитель равен произведению диагональных элементов с учётом перестановок строк: ${sign < 0 ? '−' : ''}(${Array.from({length:n}, (_,i) => fmt(a[i][i])).join(' × ')}) = ${fmt(value)}.`);
  return { value, steps: log.join('\n\n') };
}

function inverseSteps(input) {
  if (input.length !== input[0].length) throw Error('Обратная матрица существует только для квадратной матрицы.');
  const n = input.length;
  const a = input.map((row, i) => [...row, ...Array.from({length:n}, (_,j) => +(i === j))]);
  const log = [`Приписываем к A единичную матрицу:\n${showMatrix(a)}`];
  for (let i = 0; i < n; i++) {
    let p = i;
    for (let r = i + 1; r < n; r++) if (Math.abs(a[r][i]) > Math.abs(a[p][i])) p = r;
    if (Math.abs(a[p][i]) < EPS) throw Error('Матрица вырожденная: обратной матрицы не существует.');
    if (p !== i) { [a[i],a[p]] = [a[p],a[i]]; log.push(`Меняем строки ${i + 1} и ${p + 1}:\n${showMatrix(a)}`); }
    const pivot = a[i][i];
    a[i] = a[i].map(v => v / pivot);
    log.push(`Делим строку ${i + 1} на ведущий элемент ${fmt(pivot)}:\n${showMatrix(a)}`);
    for (let r = 0; r < n; r++) {
      if (r === i) continue;
      const f = a[r][i];
      if (Math.abs(f) < EPS) continue;
      a[r] = a[r].map((v,c) => v - f * a[i][c]);
      a[r] = a[r].map(v => Math.abs(v) < EPS ? 0 : v);
      log.push(`R${r + 1} ← R${r + 1} − (${fmt(f)})·R${i + 1}:\n${showMatrix(a)}`);
    }
  }
  const result = a.map(row => row.slice(n));
  log.push(`Правая часть полученной матрицы — A⁻¹:\n${showMatrix(result)}`);
  return { value: result, steps: log.join('\n\n') };
}

function rrefSteps(input) {
  const a = clone(input), rows = a.length, cols = a[0].length, log = [`Начинаем с матрицы A:\n${showMatrix(a)}`];
  let pivotRow = 0;
  for (let col = 0; col < cols && pivotRow < rows; col++) {
    let p = pivotRow;
    while (p < rows && Math.abs(a[p][col]) < EPS) p++;
    if (p === rows) { log.push(`В столбце ${col + 1} нет ведущего элемента, пропускаем его.`); continue; }
    if (p !== pivotRow) { [a[p],a[pivotRow]] = [a[pivotRow],a[p]]; log.push(`Меняем строки ${p + 1} и ${pivotRow + 1}:\n${showMatrix(a)}`); }
    const pivot = a[pivotRow][col];
    a[pivotRow] = a[pivotRow].map(v => v / pivot);
    log.push(`Делим строку ${pivotRow + 1} на ${fmt(pivot)}, чтобы ведущий элемент стал 1:\n${showMatrix(a)}`);
    for (let r = 0; r < rows; r++) {
      if (r === pivotRow) continue;
      const f = a[r][col];
      if (Math.abs(f) < EPS) continue;
      a[r] = a[r].map((v,c) => v - f * a[pivotRow][c]);
      a[r] = a[r].map(v => Math.abs(v) < EPS ? 0 : v);
      log.push(`R${r + 1} ← R${r + 1} − (${fmt(f)})·R${pivotRow + 1}:\n${showMatrix(a)}`);
    }
    pivotRow++;
  }
  return { value: a, steps: `${log.join('\n\n')}\n\nИтоговая приведённая ступенчатая матрица:\n${showMatrix(a)}` };
}

function rankSteps(input) {
  const reduced = rrefSteps(input);
  const rank = reduced.value.filter(row => row.some(v => Math.abs(v) > EPS)).length;
  return { value: rank, steps: `${reduced.steps}\n\nРанг — количество ненулевых строк после преобразований: ${rank}.` };
}

function traceSteps(a) {
  if (a.length !== a[0].length) throw Error('След вычисляется только для квадратной матрицы.');
  const diagonal = a.map((row,i) => row[i]);
  const value = diagonal.reduce((sum,x) => sum + x, 0);
  return { value, steps: `Берём элементы главной диагонали: ${diagonal.map(fmt).join(' + ')}.\nСкладываем: ${fmt(value)}.` };
}

function addSteps(a,b,subtract = false) {
  if (!sameSize(a,b)) throw Error('Для сложения и вычитания размеры матриц должны совпадать.');
  const sign = subtract ? -1 : 1;
  const operation = subtract ? 'вычитания' : 'сложения';
  const result = add(a,b,sign);
  const rows = a.map((row,i) => row.map((v,j) => `${fmt(v)} ${subtract ? '−' : '+'} ${fmt(b[i][j])} = ${fmt(result[i][j])}`).join('   ')).join('\n');
  return { value: result, steps: `Складываем/вычитаем соответствующие элементы матриц:\n${rows}\n\nРезультат ${operation}:\n${showMatrix(result)}` };
}

function multiplySteps(a,b) {
  if (a[0].length !== b.length) throw Error('Число столбцов A должно совпадать с числом строк B.');
  const result = multiply(a,b), details = [];
  for (let i = 0; i < a.length; i++) for (let j = 0; j < b[0].length; j++) {
    const terms = a[i].map((v,k) => `${fmt(v)}·${fmt(b[k][j])}`);
    details.push(`Результат[${i + 1},${j + 1}] = ${terms.join(' + ')} = ${fmt(result[i][j])}`);
  }
  return { value: result, steps: `Умножаем строку каждой матрицы A на столбец матрицы B и складываем произведения:\n${details.join('\n')}\n\nРезультат:\n${showMatrix(result)}` };
}

function transposeSteps(a) {
  const result = transpose(a);
  return { value: result, steps: `Строки исходной матрицы становятся столбцами, а столбцы — строками.\nИсходная матрица:\n${showMatrix(a)}\n\nТранспонированная матрица:\n${showMatrix(result)}` };
}

function cramersSteps(A,b) {
  const d = determinantSteps(A), D = d.value;
  const log = [`Метод Крамера для системы Ax = b.\nA:\n${showMatrix(A)}\nb = (${b.map(fmt).join(', ')})`, `\n1. Находим D = det(A):\n${d.steps}`];
  if (Math.abs(D) < EPS) {
    const rankA = rankSteps(A).value;
    const rankAug = rankSteps(A.map((row,i) => [...row,b[i]])).value;
    const answer = rankA < rankAug ? 'Система не имеет решений.' : 'Единственного решения нет: система имеет бесконечно много решений.';
    log.push(`D = 0, поэтому метод Крамера не даёт единственного решения. rank(A) = ${rankA}, rank([A|b]) = ${rankAug}. ${answer}`);
    return { answer, steps: log.join('\n\n') };
  }
  const x = [];
  for (let j = 0; j < A.length; j++) {
    const Aj = A.map((row,i) => row.map((v,c) => c === j ? b[i] : v));
    const dj = determinantSteps(Aj);
    x[j] = dj.value / D;
    log.push(`\n${j + 2}. Заменяем столбец ${j + 1} матрицы A на b:\nA${j + 1}:\n${showMatrix(Aj)}\nD${j + 1} = det(A${j + 1}) = ${fmt(dj.value)}.\nРасчёт определителя:\n${dj.steps}\nx${j + 1} = D${j + 1}/D = ${fmt(dj.value)}/${fmt(D)} = ${fmt(x[j])}.`);
  }
  log.push(`\nОтвет: x = (${x.map(fmt).join(', ')}).`);
  return { answer: `Решение: x = (${x.map(fmt).join(', ')})`, steps: log.join('\n\n') };
}

function gaussSteps(A,b) {
  if (A.length !== A[0].length) throw Error('Метод Гаусса рассчитан на квадратную матрицу A.');
  const n = A.length, m = A.map((row,i) => [...row,b[i]]), log = [`Составляем расширенную матрицу [A | b]:\n${showMatrix(m)}`];
  let pivotRow = 0;
  for (let col = 0; col < n && pivotRow < n; col++) {
    let p = pivotRow;
    for (let r = pivotRow + 1; r < n; r++) if (Math.abs(m[r][col]) > Math.abs(m[p][col])) p = r;
    if (Math.abs(m[p][col]) < EPS) { log.push(`В столбце ${col + 1} нет ведущего элемента; переходим к следующему.`); continue; }
    if (p !== pivotRow) { [m[p],m[pivotRow]] = [m[pivotRow],m[p]]; log.push(`Меняем местами строки ${pivotRow + 1} и ${p + 1}:\n${showMatrix(m)}`); }
    const pivot = m[pivotRow][col];
    log.push(`Ведущий элемент: ${fmt(pivot)} (строка ${pivotRow + 1}, столбец ${col + 1}).`);
    for (let r = pivotRow + 1; r < n; r++) {
      const f = m[r][col] / pivot;
      if (Math.abs(f) < EPS) continue;
      for (let c = col; c <= n; c++) m[r][c] -= f * m[pivotRow][c];
      m[r] = m[r].map(v => Math.abs(v) < EPS ? 0 : v);
      log.push(`Обнуляем элемент под ведущим: R${r + 1} ← R${r + 1} − (${fmt(f)})·R${pivotRow + 1}.\n${showMatrix(m)}`);
    }
    pivotRow++;
  }
  const rankA = rankSteps(A).value, rankAug = rankSteps(m.map(row => row.slice(0,n+1))).value;
  if (rankA < rankAug) {
    log.push(`rank(A) = ${rankA}, rank([A|b]) = ${rankAug}. Ранги различны — система несовместна.`);
    return { answer: 'Система не имеет решений.', steps: log.join('\n\n') };
  }
  if (rankA < n) {
    log.push(`rank(A) = rank([A|b]) = ${rankA} < ${n}. Есть свободные переменные, решений бесконечно много.`);
    return { answer: 'Система имеет бесконечно много решений.', steps: log.join('\n\n') };
  }
  const x = Array(n).fill(0);
  log.push('Обратный ход: подставляем найденные значения снизу вверх.');
  for (let r = n - 1; r >= 0; r--) {
    let sum = m[r][n];
    for (let c = r + 1; c < n; c++) sum -= m[r][c] * x[c];
    x[r] = sum / m[r][r];
    log.push(`x${r + 1} = ${fmt(x[r])}`);
  }
  log.push(`Ответ: x = (${x.map(fmt).join(', ')}).`);
  return { answer: `Решение: x = (${x.map(fmt).join(', ')})`, steps: log.join('\n\n') };
}

function run(event) {
  try {
    const op = event.currentTarget.dataset.op;
    let result, answer;
    if (op === 'cramer' || op === 'gauss') {
      const {A,b} = parseSystem();
      result = op === 'cramer' ? cramersSteps(A,b) : gaussSteps(A,b);
      output.textContent = result.answer;
      stepsOutput.textContent = result.steps;
      return;
    }
    const a = parse($('matrixA').value), b = parse($('matrixB').value);
    switch (op) {
      case 'add': result=addSteps(a,b); answer=showMatrix(result.value); break;
      case 'subtract': result=addSteps(a,b,true); answer=showMatrix(result.value); break;
      case 'multiply': result=multiplySteps(a,b); answer=showMatrix(result.value); break;
      case 'transposeA': result=transposeSteps(a); answer=showMatrix(result.value); break;
      case 'transposeB': result=transposeSteps(b); answer=showMatrix(result.value); break;
      case 'detA': result=determinantSteps(a); answer=`Определитель A = ${fmt(result.value)}`; break;
      case 'detB': result=determinantSteps(b); answer=`Определитель B = ${fmt(result.value)}`; break;
      case 'inverseA': result=inverseSteps(a); answer=showMatrix(result.value); break;
      case 'inverseB': result=inverseSteps(b); answer=showMatrix(result.value); break;
      case 'rankA': result=rankSteps(a); answer=`Ранг A = ${result.value}`; break;
      case 'traceA': result=traceSteps(a); answer=`След A = ${fmt(result.value)}`; break;
      case 'rrefA': result=rrefSteps(a); answer=showMatrix(result.value); break;
      default: throw Error('Неизвестная операция.');
    }
    output.textContent = answer;
    stepsOutput.textContent = result.steps;
  } catch (error) {
    output.textContent = error.message || 'Ошибка вычисления.';
    stepsOutput.textContent = '';
  }
}

document.querySelectorAll('[data-op]').forEach(button => button.addEventListener('click',run));
$('sample').addEventListener('click',() => {
  $('matrixA').value='2 1 0\n1 3 1\n0 1 2';
  $('matrixB').value='1\n2\n3';
  output.textContent='Пример загружен. Выберите действие.';
  stepsOutput.textContent='';
});
