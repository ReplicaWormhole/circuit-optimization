#!/usr/bin/env python3
"""Render complete post-gate matrices of the accepted four-qubit circuit."""
import argparse
import colorsys
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from check_circuit import cnot, embedded_one_qubit, evaluate, one_qubit_matrix, right_shift

BG, FG, MUTED, GOLD = (12,18,30), (230,238,248), (147,166,189), (255,201,89)

@lru_cache(maxsize=None)
def font(size, mono=False):
    name = 'DejaVuSansMono.ttf' if mono else 'DejaVuSans.ttf'
    return ImageFont.truetype('/usr/share/fonts/truetype/dejavu/' + name, size)


def prefixes(candidate):
    evaluate(candidate)  # Shared input validation.
    n = candidate['n']
    u = np.eye(1 << n, dtype=complex)
    units = [u.copy()]
    for gate in candidate['gates']:
        g = (cnot(gate['control'], gate['target'], n) if gate['gate'].lower() == 'cx'
             else embedded_one_qubit(one_qubit_matrix(gate), gate['qubit'], n))
        u = g @ u
        units.append(u.copy())
    v = right_shift(n)
    return np.array(units), np.array([u @ v @ u.conj().T for u in units])


def entry_text(z):
    if abs(z) < 1e-10:
        return '0'
    if abs(z.imag) < 1e-10:
        return f'{z.real:.3f}'
    if abs(z.real) < 1e-10:
        return f'{z.imag:.3f}i'
    return f'{z.real:.3f}\n{z.imag:+.3f}i'


def color(z):
    magnitude = min(1., abs(z))
    if magnitude < 1e-10:
        return (18,26,41)
    rgb = colorsys.hsv_to_rgb((np.angle(z)/(2*np.pi)) % 1, .55, .18+.39*magnitude)
    return tuple(round(255*c) for c in rgb)


def gate_label(gate):
    name = gate['gate'].upper()
    if name == 'CX':
        return f"CX q{gate['control']} -> q{gate['target']}"
    keys = ('theta','phi','lam') if name == 'U3' else ('theta',)
    args = ', '.join(str(gate[k]) for k in keys if k in gate)
    return f"{name}({args}) on q{gate['qubit']}" if args else f"{name} on q{gate['qubit']}"


def draw_matrix(draw, matrix, x, title):
    y, cw, ch = 650, 96, 80
    draw.text((x,540),title,font=font(36),fill=FG)
    draw.text((x,590),'columns: input basis     rows: output basis',font=font(23),fill=MUTED)
    for b in range(16):
        draw.text((x+cw*(b+.5),y-28),f'{b:04b}',font=font(23,True),fill=MUTED,anchor='mm')
        draw.text((x-16,y+ch*(b+.5)),f'{b:04b}',font=font(23,True),fill=MUTED,anchor='rm')
    for row in range(16):
        for col in range(16):
            z = matrix[row,col]
            draw.rectangle((x+col*cw,y+row*ch,x+(col+1)*cw-2,y+(row+1)*ch-2),fill=color(z))
            draw.multiline_text((x+(col+.5)*cw,y+(row+.5)*ch),entry_text(z),font=font(21,True),
                                fill=FG if abs(z)>1e-10 else (88,104,126),anchor='mm',align='center',spacing=4)


def draw_circuit(draw, gates, step):
    for q in range(4):
        y = 235+70*q
        draw.text((100,y),f'q{q}',font=font(28,True),fill=FG,anchor='mm')
        draw.line((160,y,3710,y),fill=(69,85,108),width=3)
    for k,gate in enumerate(gates,1):
        x = 220+(k-1)*116
        ink = GOLD if k == step else ((96,210,200) if k < step else MUTED)
        if k == step:
            draw.rounded_rectangle((x-52,170,x+52,492),radius=12,outline=GOLD,width=4)
        draw.text((x,186),str(k),font=font(21,True),fill=ink,anchor='mm')
        if gate['gate'].lower() == 'cx':
            yc,yt = 235+70*gate['control'],235+70*gate['target']
            draw.line((x,yc,x,yt),fill=ink,width=4)
            draw.ellipse((x-9,yc-9,x+9,yc+9),fill=ink)
            draw.ellipse((x-18,yt-18,x+18,yt+18),fill=BG,outline=ink,width=3)
            draw.line((x-12,yt,x+12,yt),fill=ink,width=3)
            draw.line((x,yt-12,x,yt+12),fill=ink,width=3)
        else:
            y = 235+70*gate['qubit']
            draw.rounded_rectangle((x-43,y-25,x+43,y+25),radius=6,fill=BG,outline=ink,width=3)
            draw.text((x,y),gate['gate'].upper(),font=font(24),fill=ink,anchor='mm')
            if 'theta' in gate:
                angle = str(gate['theta'])
                if gate['gate'].lower() == 'u3':
                    angle = 'λ='+str(gate['lam'])
                draw.text((x,y+36),angle.replace('*pi','π').replace('pi','π'),font=font(19),fill=ink,anchor='mm')
    draw.text((220,506),'Chronological gates, left to right. Gold = gate just applied. Full U3 parameters appear in the header.',font=font(23),fill=MUTED)


def render(candidate, units, shifted, step):
    image = Image.new('RGB',(3840,2160),BG)
    draw = ImageDraw.Draw(image)
    draw.text((100,50),'FOUR-QUBIT SHIFT  /  11-CNOT DIAGONALIZER',font=font(45),fill=FG)
    label = 'Initial state: U₀ = I' if step == 0 else gate_label(candidate['gates'][step-1])
    draw.text((100,116),f'Step {step:02d} / {len(candidate["gates"])}     {label}',font=font(34),fill=GOLD)
    draw_circuit(draw,candidate['gates'],step)
    draw_matrix(draw,units[step],220,'Cumulative circuit Uₖ = Gₖ ... G₂ G₁')
    draw_matrix(draw,shifted[step],2160,'Shift in the evolving basis: Uₖ V₄ Uₖ†')
    off = shifted[step]-np.diag(np.diag(shifted[step]))
    draw.text((2160,1960),f'Maximum |off-diagonal entry| = {np.max(np.abs(off)):.3e}',font=font(28),fill=GOLD)
    draw.text((220,1960),'Full 16 × 16 matrices; binary basis 0000 ... 1111; q0 is the most significant bit.',font=font(24),fill=MUTED)
    draw.text((220,2010),'Entries rounded to 3 decimals. |z| < 1e-10 displays as zero. All global phases retained.',font=font(24),fill=MUTED)
    draw.text((220,2060),'Color: phase (hue) + magnitude (brightness).',font=font(25),fill=MUTED)
    for j,(label,z) in enumerate([('+1',1),('+i',1j),('−1',-1),('−i',-1j)]):
        x = 2160+j*240
        draw.rectangle((x,2036,x+50,2086),fill=color(z))
        draw.text((x+65,2061),label,font=font(28),fill=FG,anchor='lm')
    return image


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=ROOT/'experiments/runs/411/candidate.json')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'video')
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text())
    if candidate['n'] != 4 or len(candidate['gates']) > 30:
        parser.error('This 4K layout supports four qubits and at most 30 gates.')
    if args.output.exists():
        parser.error('Output exists; choose a fresh --output to preserve existing renders.')
    units,shifted = prefixes(candidate)
    args.output.mkdir(parents=True)
    frames = args.output/'frames'
    frames.mkdir()
    np.savez_compressed(args.output/'matrices.npz',U=units,conjugated_shift=shifted)
    for step in range(len(units)):
        render(candidate,units,shifted,step).save(frames/f'step_{step:03d}.png')
        print(f'Rendered step {step}/{len(units)-1}',flush=True)
    report = {'candidate':str(args.candidate.resolve()),'candidate_sha256':hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
              'gate_count':len(candidate['gates']),'matrix_states':len(units),'resolution':[3840,2160],
              'seconds_per_state':2,'extra_final_hold_seconds':3,'checker':evaluate(candidate),
              'render_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'display_scope':'Floating-point visualization, 3 decimal labels, threshold 1e-10; no new exact proof.',
              'numpy_version':np.__version__}
    (args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__ == '__main__':
    main()
