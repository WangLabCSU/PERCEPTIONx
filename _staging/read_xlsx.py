import zipfile, re, sys, csv, os
from xml.etree import ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
path = r"C:\Users\Lenovo\Desktop\PERCEPTION\_staging\Maynard_Supp2_Demo.xlsx"
outdir = r"C:\Users\Lenovo\Desktop\PERCEPTION\_staging"

z = zipfile.ZipFile(path)
# shared strings
shared = []
if 'xl/sharedStrings.xml' in z.namelist():
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    for si in root.findall('m:si', NS):
        txt = ''.join(t.text or '' for t in si.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t'))
        shared.append(txt)

# sheet name -> rId mapping
wb = ET.fromstring(z.read('xl/workbook.xml'))
rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
NS_R = {'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
sheet_rid = {}
for sh in wb.findall('m:sheets/m:sheet', NS):
    sheet_rid[sh.get('name')] = sh.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
rid_path = {}
for rel in rels.findall('{http://schemas.openxmlformats.org/package/2006/relationships}Relationship'):
    rid_path[rel.get('Id')] = rel.get('Target')

def colnum(ref):
    m = re.match(r'([A-Z]+)(\d+)', ref)
    col = 0
    for ch in m.group(1):
        col = col * 26 + (ord(ch) - 64)
    return col, int(m.group(2))

for sname, rid in sheet_rid.items():
    tgt = rid_path[rid]
    if not tgt.startswith('xl/'):
        tgt = 'xl/' + tgt
    root = ET.fromstring(z.read(tgt))
    rows = []
    for row in root.iter('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row'):
        cells = {}
        maxc = 0
        for c in row.findall('m:c', NS):
            ref = c.get('r')
            if not ref: continue
            ci, _ = colnum(ref)
            t = c.get('t')
            v = c.find('m:v', NS)
            val = ''
            if v is not None and v.text is not None:
                if t == 's':
                    val = shared[int(v.text)]
                else:
                    val = v.text
            cells[ci] = val
            maxc = max(maxc, ci)
        rows.append([cells.get(i, '') for i in range(1, maxc + 1)])
    fname = os.path.join(outdir, 'sheet_' + re.sub(r'\W+', '_', sname) + '.csv')
    with open(fname, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerows(rows)
    print(f"sheet '{sname}' -> {fname} ({len(rows)} rows)")
print("DONE")
