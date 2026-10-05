#!/usr/bin/env python3
"""tab_edit.py - cau_tab 表点编辑器（简版）：只改值，z 网格不动。

用法: python tab_edit.py <in_table> <out_table> SPEC [SPEC ...]
SPEC:  sheet/site/op/zlo/zhi/amp      (site 可为 top|bridge|hollow|all)
  op=bump  : E'(z) = E(z) + amp * sin^2(pi*(z-zlo)/(zhi-zlo))   C1 钟形（位垒下压/阱上抬/通道平移）
  op=ramp0 : E'(z) = E(z) * (1 - smoothstep(z))                 远场平滑归零（amp 忽略；用于 E0）
  op=flat  : E'(z) = 平滑趋近窗口末点值（去远场梯度；amp 忽略；适用于本身渐平的表）
  op=hold  : E'(z) = 平滑锁定到窗口起点值（远场保位——4 A 外孤立原子，通道应无 z 依赖；
             current 谱系 h/U 远场持续下降属导出伪影，用本操作截断）
窗口外点不动。输出表顶部追加 '# edit:' 注记行（随装载回显进运行记录）。
"""
import sys, math

def load(path):
    lines, sheet, site = [], None, None
    for raw in open(path, encoding='utf-8'):
        s = raw.rstrip('\n').rstrip('\r')
        p = s.split()
        if not p or p[0].startswith('#'):
            lines.append(('c', s, None, None)); continue
        if p[0] == 'LATTICE':
            lines.append(('c', s, None, None)); continue
        if p[0] == 'SHEET':
            sheet = p[1]; lines.append(('h', s, None, None)); continue
        if p[0] == 'SITE':
            site = p[1]; lines.append(('h', s, None, None)); continue
        try:
            z, e = float(p[0]), float(p[1])
            lines.append(('d', s, sheet, site, z, e))
        except ValueError:
            lines.append(('c', s, None, None))
    return lines

def bell(z, zlo, zhi):
    t = min(max((z - zlo) / (zhi - zlo), 0.0), 1.0)
    return math.sin(math.pi * t) ** 2

def smoothstep(z, zlo, zhi):
    t = min(max((z - zlo) / (zhi - zlo), 0.0), 1.0)
    return 3 * t * t - 2 * t * t * t

def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    specs = []
    for a in sys.argv[3:]:
        sh, site, op, zlo, zhi, amp = a.split('/')
        specs.append((sh, site, op, float(zlo), float(zhi), float(amp)))
    lines = load(src)
    # flat/hold 参考值：flat 取窗口内最大 z 点的值，hold 取窗口内最小 z 点的值
    ref = {}
    for sh, st, op, zlo, zhi, amp in specs:
        if op not in ('flat', 'hold'):
            continue
        pick_max = (op == 'flat')
        for ln in lines:
            if ln[0] == 'd':
                _, s, sheet, site, z, e = ln
                if sh == sheet and (st == site or st == 'all') and zlo <= z <= zhi:
                    k = (sheet, site)
                    if k not in ref or (z > ref[k][0]) == pick_max:
                        ref[k] = (z, e)
    out, changed = [], []
    hit = [0] * len(specs)
    for ln in lines:
        if ln[0] == 'd':
            _, s, sheet, site, z, e = ln
            e2 = e
            for k, (sh, st, op, zlo, zhi, amp) in enumerate(specs):
                if sh != sheet or not (st == site or st == 'all'):
                    continue
                if not (zlo <= z <= zhi):
                    continue
                hit[k] += 1
                if op == 'bump':
                    e2 += amp * bell(z, zlo, zhi)
                elif op == 'ramp0':
                    e2 *= (1.0 - smoothstep(z, zlo, zhi))
                elif op == 'flat':
                    e2 += (ref[(sheet, site)][1] - e) * smoothstep(z, zlo, zhi)
                elif op == 'hold':
                    e2 += (ref[(sheet, site)][1] - e) * smoothstep(z, zlo, zhi)
                else:
                    sys.exit('unknown op: ' + op)
            if e2 != e:
                changed.append('%s/%s z=%.2f %+.3f->%+.3f' % (sheet, site, z, e, e2))
            out.append('  %9.6f  %24.16e' % (z, e2))
        else:
            out.append(ln[1])
    with open(dst, 'w', encoding='utf-8', newline='\n') as f:
        f.write('# edits (tab_edit.py, base=%s):\n' % src)
        for a in sys.argv[3:]:
            f.write('# edit: %s\n' % a)
        f.write('\n'.join(out) + '\n')
    bad = [('/'.join(map(str, specs[k])) , n) for k, n in enumerate(hit) if n == 0]
    if bad:
        for s, n in bad:
            sys.stderr.write('SPEC matched no points (sheet/site name or window wrong): %s\n' % s)
        sys.exit(1)
    print('points changed: %d' % len(changed))
    for c in changed:
        print('  ' + c)

if __name__ == '__main__':
    main()
