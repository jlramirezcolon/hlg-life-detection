"""Extract one face of a TrueType/OpenType collection (.ttc) into a standalone font file, without fontTools.

A .ttc holds several fonts whose table directories point into shared table data; copying one font's tables into
a new file with a fresh table directory gives a standalone font that matplotlib/FreeType can load.
Usage: python3 extract_ttc_face.py <collection.ttc> "<full font name>" <out.ttf>
       python3 extract_ttc_face.py <collection.ttc>            (lists the faces)
"""
import struct
import sys


def faces(data):
    assert data[:4] == b'ttcf', 'not a font collection'
    n = struct.unpack('>I', data[8:12])[0]
    return list(struct.unpack('>%dI' % n, data[12:12 + 4 * n]))


def tables(data, off):
    num = struct.unpack('>H', data[off + 4:off + 6])[0]
    recs = []
    for i in range(num):
        tag, chk, o, ln = struct.unpack('>4sIII', data[off + 12 + 16 * i:off + 28 + 16 * i])
        recs.append((tag, chk, o, ln))
    return data[off:off + 4], recs


def full_name(data, off):
    _, recs = tables(data, off)
    t = next(r for r in recs if r[0] == b'name')
    base = t[2]
    fmt, count, soff = struct.unpack('>HHH', data[base:base + 6])
    best = None
    for i in range(count):
        pid, eid, lid, nid, ln, o = struct.unpack('>6H', data[base + 6 + 12 * i:base + 18 + 12 * i])
        if nid != 4:
            continue
        raw = data[base + soff + o:base + soff + o + ln]
        if pid == 3 and lid == 0x409:
            return raw.decode('utf-16-be')
        if pid == 1 and best is None:
            best = raw.decode('mac_roman')
    return best


def extract(data, off):
    version, recs = tables(data, off)
    n = len(recs)
    es = max(k for k in range(16) if 2 ** k <= n)
    sr = 16 * 2 ** es
    head = version + struct.pack('>HHHH', n, sr, es, 16 * n - sr)
    pos = 12 + 16 * n
    dir_, body = b'', b''
    for tag, chk, o, ln in sorted(recs):
        blob = data[o:o + ln]
        dir_ += struct.pack('>4sIII', tag, chk, pos + len(body), ln)
        body += blob + b'\0' * ((4 - ln % 4) % 4)
    return head + dir_ + body


if __name__ == '__main__':
    data = open(sys.argv[1], 'rb').read()
    offs = faces(data)
    if len(sys.argv) == 2:
        for i, o in enumerate(offs):
            print(i, full_name(data, o))
        sys.exit()
    want = sys.argv[2]
    off = next((o for o in offs if full_name(data, o) == want), None)
    assert off is not None, f'{want} not found'
    open(sys.argv[3], 'wb').write(extract(data, off))
    print('wrote', sys.argv[3])
