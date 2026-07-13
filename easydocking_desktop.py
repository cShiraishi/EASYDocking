"""
EASYDocking Desktop — Molecular Docking Platform
PyQt6 local application — mirrors the full Streamlit web app.
"""
# ── rdkit.six compatibility patch (meeko 0.5.0) ────────────────────────────
import sys, io, types
try:
    from rdkit import six
except ImportError:
    _six = types.ModuleType('rdkit.six')
    _six.StringIO = io.StringIO
    _six.BytesIO  = io.BytesIO
    sys.modules['rdkit.six'] = _six

import os, shutil, subprocess, tempfile, zipfile, traceback, warnings
warnings.filterwarnings('ignore')


def resolve_babel_bin(name):
    """Resolve o caminho de um executável do OpenBabel (obabel/obrms).

    Procura primeiro no PATH; se não achar (comum quando o app roda sob o
    Anaconda/GUI sem o /opt/homebrew/bin no PATH), tenta locais de instalação
    conhecidos (Homebrew Apple Silicon/Intel, conda, Linux). Retorna o nome
    original como último recurso para preservar a mensagem de erro padrão.
    """
    found = shutil.which(name)
    if found:
        return found
    candidates = [
        f"/opt/homebrew/bin/{name}",          # Homebrew (Apple Silicon)
        f"/usr/local/bin/{name}",             # Homebrew (Intel) / Linux
        os.path.join(os.path.dirname(sys.executable), name),  # conda env
        f"/usr/bin/{name}",
    ]
    for path in candidates:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return name


# Executáveis do OpenBabel resolvidos uma vez na inicialização
OBABEL = resolve_babel_bin("obabel")
OBRMS = resolve_babel_bin("obrms")

import numpy as np
import pandas as pd
import requests

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QPushButton, QFileDialog, QTextEdit, QProgressBar,
    QLineEdit, QGroupBox, QRadioButton, QTableWidget, QTableWidgetItem,
    QHeaderView, QDoubleSpinBox, QCheckBox, QComboBox, QSpinBox,
    QFrame, QSplitter, QScrollArea, QMessageBox, QStatusBar,
    QButtonGroup, QAbstractItemView, QSizePolicy,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QUrl
from PyQt6.QtGui  import QFont, QColor, QPixmap
from PyQt6.QtWebEngineWidgets import QWebEngineView

import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# ── Palette ────────────────────────────────────────────────────────────────────
C = {
    'bg':       '#0b0f14', 'surface':  '#111720', 'surface2': '#161d28',
    'surface3': '#1c2535', 'border':   '#1e2d3d', 'border2':  '#243044',
    'accent':   '#3d8ef0', 'accent_h': '#5aa3f7', 'accent_d': '#2a6bcc',
    'purple':   '#7c3aed', 'green':    '#22c55e', 'amber':    '#f59e0b',
    'red':      '#ef4444', 'text':     '#d4dce8', 'text2':    '#8494a8',
    'text3':    '#4a5568',
}

QSS = f"""
QMainWindow,QDialog,QWidget{{background:{C['bg']};color:{C['text']};
    font-family:"Inter","Segoe UI",Arial,sans-serif;font-size:12px;}}
QTabWidget::pane{{border:1px solid {C['border']};background:{C['surface']};
    border-bottom-left-radius:8px;border-bottom-right-radius:8px;}}
QTabBar{{background:transparent;}}
QTabBar::tab{{background:transparent;color:{C['text2']};padding:10px 22px;
    font-size:12px;font-weight:500;border:none;
    border-bottom:2px solid transparent;margin-right:2px;}}
QTabBar::tab:hover{{color:{C['text']};background:{C['surface2']};
    border-bottom:2px solid {C['border2']};}}
QTabBar::tab:selected{{color:{C['accent']};background:{C['surface']};
    border-bottom:2px solid {C['accent']};font-weight:600;}}
QPushButton{{background:{C['surface3']};color:{C['text']};
    border:1px solid {C['border2']};border-radius:6px;
    padding:6px 14px;font-size:12px;font-weight:500;}}
QPushButton:hover{{background:{C['surface2']};border-color:{C['accent']};color:#f0f4f9;}}
QPushButton:pressed{{background:{C['bg']};}}
QPushButton:disabled{{background:{C['surface']};color:{C['text3']};
    border-color:{C['border']};}}
QLineEdit,QTextEdit{{background:{C['surface2']};color:{C['text']};
    border:1px solid {C['border2']};border-radius:6px;padding:6px 10px;}}
QLineEdit:focus,QTextEdit:focus{{border-color:{C['accent']};}}
QDoubleSpinBox,QSpinBox{{background:{C['surface2']};color:{C['text']};
    border:1px solid {C['border2']};border-radius:6px;padding:4px 8px;}}
QDoubleSpinBox:focus,QSpinBox:focus{{border-color:{C['accent']};}}
QComboBox{{background:{C['surface2']};color:{C['text']};
    border:1px solid {C['border2']};border-radius:6px;
    padding:5px 10px;min-height:28px;}}
QComboBox:hover,QComboBox:focus{{border-color:{C['accent']};}}
QComboBox QAbstractItemView{{background:{C['surface2']};color:{C['text']};
    border:1px solid {C['border2']};
    selection-background-color:{C['accent_d']};}}
QGroupBox{{background:{C['surface']};border:1px solid {C['border']};
    border-radius:8px;margin-top:14px;padding-top:8px;
    font-size:10px;font-weight:700;color:{C['text2']};
    letter-spacing:0.5px;}}
QGroupBox::title{{subcontrol-origin:margin;subcontrol-position:top left;
    padding:0 8px;left:12px;background:{C['bg']};
    color:{C['accent']};font-size:10px;font-weight:700;letter-spacing:1px;}}
QTableWidget{{background:{C['surface']};color:{C['text']};
    border:1px solid {C['border']};border-radius:8px;
    gridline-color:{C['border']};
    alternate-background-color:{C['surface2']};font-size:11px;}}
QTableWidget::item{{padding:5px 8px;border:none;}}
QTableWidget::item:selected{{background:{C['accent_d']};color:#f0f4f9;}}
QTableWidget::item:hover{{background:{C['surface3']};}}
QHeaderView::section{{background:{C['surface2']};color:{C['text2']};
    border:none;border-right:1px solid {C['border']};
    border-bottom:1px solid {C['border']};
    padding:6px 10px;font-size:10px;font-weight:700;
    letter-spacing:0.5px;}}
QProgressBar{{background:{C['surface2']};border:1px solid {C['border']};
    border-radius:4px;text-align:center;
    color:{C['text2']};font-size:10px;}}
QProgressBar::chunk{{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,
    stop:0 {C['accent_d']},stop:1 {C['accent']});border-radius:4px;}}
QScrollBar:vertical{{background:{C['surface']};width:8px;border-radius:4px;margin:0;}}
QScrollBar::handle:vertical{{background:{C['border2']};border-radius:4px;min-height:24px;}}
QScrollBar::add-line:vertical,QScrollBar::sub-line:vertical{{height:0;}}
QScrollBar:horizontal{{background:{C['surface']};height:8px;border-radius:4px;}}
QScrollBar::handle:horizontal{{background:{C['border2']};border-radius:4px;min-width:24px;}}
QScrollBar::add-line:horizontal,QScrollBar::sub-line:horizontal{{width:0;}}
QCheckBox,QRadioButton{{color:{C['text']};spacing:8px;font-size:12px;}}
QCheckBox::indicator,QRadioButton::indicator{{width:16px;height:16px;
    border:1px solid {C['border2']};border-radius:4px;background:{C['surface2']};}}
QCheckBox::indicator:checked{{background:{C['accent']};border-color:{C['accent']};}}
QRadioButton::indicator{{border-radius:8px;}}
QRadioButton::indicator:checked{{background:{C['accent']};border-color:{C['accent']};}}
QStatusBar{{background:{C['surface']};color:{C['text2']};
    border-top:1px solid {C['border']};font-size:11px;padding:2px 8px;}}
QScrollArea{{border:none;background:transparent;}}
QSplitter::handle{{background:{C['border']};width:1px;height:1px;}}
QToolTip{{background:{C['surface3']};color:{C['text']};
    border:1px solid {C['border2']};border-radius:6px;padding:5px 9px;font-size:11px;}}
"""

# ── Console widget ─────────────────────────────────────────────────────────────
class Console(QTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True); self.setMaximumHeight(150); self.setMinimumHeight(70)
        self.setStyleSheet(f"""QTextEdit{{background:{C['bg']};color:{C['accent']};
            font-family:'JetBrains Mono','Fira Code','Consolas',monospace;
            font-size:11px;border:1px solid {C['border']};border-radius:6px;padding:6px;}}""")
    def _ts(self):
        from datetime import datetime; return datetime.now().strftime('%H:%M:%S')
    def info(self, m):    self.append(f"<span style='color:{C['text3']}'>{self._ts()}</span> <span style='color:{C['text2']}'>{m}</span>")
    def ok(self, m):      self.append(f"<span style='color:{C['text3']}'>{self._ts()}</span> <span style='color:{C['green']}'>✔ {m}</span>")
    def warn(self, m):    self.append(f"<span style='color:{C['text3']}'>{self._ts()}</span> <span style='color:{C['amber']}'>⚠ {m}</span>")
    def err(self, m):     self.append(f"<span style='color:{C['text3']}'>{self._ts()}</span> <span style='color:{C['red']}'>✖ {m}</span>")


# ── 3D Viewer (py3Dmol via QWebEngineView) ─────────────────────────────────────
class MolViewer(QWebEngineView):
    """Embeds py3Dmol interactive viewer via generated HTML."""
    _EMPTY_HTML = f"""<!DOCTYPE html><html><body
        style="margin:0;background:{C['bg']};display:flex;align-items:center;
               justify-content:center;height:100vh;color:{C['text2']};
               font-family:sans-serif;font-size:13px;">
        <span>Load a structure to visualize</span></body></html>"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(400)
        self.setHtml(self._EMPTY_HTML)

    def show_structure(self, pdbqt_or_pdb: str, fmt: str = 'pdbqt',
                       ligand_pdbqt: str = None,
                       center: tuple = None, box_size: tuple = None,
                       ref_pdbqt: str = None):
        """Render protein ± ligand ± grid box using py3Dmol JS."""
        import json
        receptor_js = json.dumps(pdbqt_or_pdb)
        html_parts = ["""<!DOCTYPE html><html><head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>body{margin:0;background:#0b0f14;}#viewer{width:100vw;height:100vh;}</style>
</head><body><div id="viewer"></div><script>
let viewer = $3Dmol.createViewer("viewer",{backgroundColor:"0x0b0f14"});
"""]
        html_parts.append(f"viewer.addModel({receptor_js},'{fmt}');\n")
        html_parts.append("viewer.setStyle({},{cartoon:{color:'spectrum',opacity:0.85}});\n")

        if ligand_pdbqt:
            lig_js = json.dumps(ligand_pdbqt)
            html_parts.append(f"viewer.addModel({lig_js},'pdbqt');\n")
            html_parts.append("viewer.setStyle({model:-1},{stick:{colorscheme:'greenCarbon'}});\n")

        if ref_pdbqt:
            ref_js = json.dumps(ref_pdbqt)
            html_parts.append(f"viewer.addModel({ref_js},'pdbqt');\n")
            html_parts.append("viewer.setStyle({model:-1},{stick:{colorscheme:'cyanCarbon',opacity:0.7}});\n")

        if center and box_size:
            cx, cy, cz = center
            sx, sy, sz = box_size
            html_parts.append(
                f"viewer.addBox({{center:{{x:{cx},y:{cy},z:{cz}}},"
                f"dimensions:{{w:{sx},h:{sy},d:{sz}}},"
                f"color:'0x22c55e',opacity:0.25}});\n"
                f"viewer.addSphere({{center:{{x:{cx},y:{cy},z:{cz}}},"
                f"radius:0.8,color:'0xef4444',opacity:0.9}});\n"
            )

        html_parts.append("viewer.zoomTo();viewer.render();\n</script></body></html>")
        self.setHtml(''.join(html_parts))

    def clear(self):
        self.setHtml(self._EMPTY_HTML)


# ── Docking Worker ─────────────────────────────────────────────────────────────
class DockingWorker(QThread):
    progress  = pyqtSignal(int, int, str)   # (done, total, ligand_name)
    result    = pyqtSignal(object)           # one result dict per ligand
    finished  = pyqtSignal()
    error     = pyqtSignal(str)

    def __init__(self, protein_path: str, ligands: list, params: dict):
        super().__init__()
        self.protein_path = protein_path
        self.ligands      = ligands   # list of (name, mol)
        self.params       = params

    def run(self):
        try:
            from vina import Vina
            from rdkit.Chem import AllChem
            from rdkit.Geometry import Point3D
            from meeko import MoleculePreparation

            p = self.params
            v = Vina(sf_name='vina')
            v.set_receptor(self.protein_path)
            v.compute_vina_maps(
                center=[p['cx'], p['cy'], p['cz']],
                box_size=[p['sx'], p['sy'], p['sz']]
            )

            total = len(self.ligands)
            for idx, (name, mol) in enumerate(self.ligands):
                self.progress.emit(idx, total, name)

                # 3D conformer
                if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                    try: AllChem.UFFOptimizeMolecule(mol)
                    except: pass
                    # Center on grid
                    conf  = mol.GetConformer()
                    coords = conf.GetPositions()
                    centroid = coords.mean(axis=0)
                    dx, dy, dz = p['cx']-centroid[0], p['cy']-centroid[1], p['cz']-centroid[2]
                    for k in range(mol.GetNumAtoms()):
                        pos = conf.GetAtomPosition(k)
                        conf.SetAtomPosition(k, Point3D(pos.x+dx, pos.y+dy, pos.z+dz))
                elif p.get('minimize'):
                    try: AllChem.MMFFOptimizeMolecule(mol)
                    except:
                        try: AllChem.UFFOptimizeMolecule(mol)
                        except: pass

                # Prepare PDBQT
                prep = MoleculePreparation()
                prep.prepare(mol)
                lig_pdbqt = prep.write_pdbqt_string()
                if not lig_pdbqt:
                    self.result.emit({'Ligand': name, 'error': 'Meeko prep failed'})
                    continue

                # Bounds check
                conf   = mol.GetConformer()
                coords = conf.GetPositions()
                c_min, c_max = coords.min(axis=0), coords.max(axis=0)
                b_min = np.array([p['cx']-p['sx']/2, p['cy']-p['sy']/2, p['cz']-p['sz']/2])
                b_max = np.array([p['cx']+p['sx']/2, p['cy']+p['sy']/2, p['cz']+p['sz']/2])
                if np.any(c_min < b_min - 15.0) or np.any(c_max > b_max + 15.0):
                    self.result.emit({'Ligand': name,
                                      'error': 'Ligand outside grid bounds (>15 Å)'})
                    continue

                # Dock
                v.set_ligand_from_string(lig_pdbqt)
                v.dock(exhaustiveness=p['exhaustiveness'], n_poses=p['n_poses'])
                energies = v.energies(n_poses=1)
                affinity = float(energies[0][0]) if len(energies) > 0 else 0.0

                t = tempfile.NamedTemporaryFile(delete=False, suffix='.pdbqt')
                t.close()
                v.write_poses(t.name, n_poses=p['n_poses'], overwrite=True)
                with open(t.name) as fh: poses_pdbqt = fh.read()
                os.unlink(t.name)

                # RMSD for redocking
                rmsd = '-'
                if '_Redocking' in name:
                    try:
                        with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as tr:
                            from rdkit import Chem
                            w = Chem.SDWriter(tr.name); w.write(mol); w.close()
                            ref_path = tr.name
                        with tempfile.NamedTemporaryFile(suffix='.pdbqt', mode='w', delete=False) as tp:
                            tp.write(poses_pdbqt); tp_path = tp.name
                        r = subprocess.run([OBRMS, '-f', ref_path, tp_path],
                                           capture_output=True, text=True)
                        for line in r.stdout.splitlines():
                            if 'RMSD' in line:
                                rmsd = round(float(line.split()[-1]), 3); break
                        os.unlink(ref_path); os.unlink(tp_path)
                    except: pass

                heavy = mol.GetNumHeavyAtoms()
                self.result.emit({
                    'Ligand':               name,
                    'Affinity (kcal/mol)':  round(affinity, 3),
                    'RMSD (Å)':             rmsd,
                    'Heavy Atoms':          heavy,
                    'Ligand Efficiency':    round(abs(affinity)/heavy, 3) if heavy else 0,
                    'PDBQT':                poses_pdbqt,
                })
            self.finished.emit()
        except Exception:
            self.error.emit(traceback.format_exc())


# ── Grid auto-calculation Worker ───────────────────────────────────────────────
class GridWorker(QThread):
    finished = pyqtSignal(float, float, float, float, float, float, str)
    error    = pyqtSignal(str)

    def __init__(self, pdb_id: str):
        super().__init__()
        self.pdb_id = pdb_id.strip().upper()

    def run(self):
        try:
            IGNORE = {"HOH","DOD","WAT","NA","CL","K","MG","CA","ZN","CU","FE",
                      "MN","CO","NI","I","BR","SO4","PO4","NO3","CO3","ACT","FMT",
                      "ACE","GOL","PEG","EDO","DMS","PG4","PGE","PE4","BME","DTT"}
            # 1. Find ligand comp_id via GraphQL
            q = '{ entry(entry_id:"%s"){ rcsb_binding_affinity{comp_id} nonpolymer_entities{nonpolymer_comp{chem_comp{id}}} } }' % self.pdb_id
            resp = requests.post('https://data.rcsb.org/graphql', json={'query': q}, timeout=15)
            comp_id = None
            if resp.ok:
                entry = resp.json().get('data', {}).get('entry', {})
                affs = entry.get('rcsb_binding_affinity')
                if affs:
                    comp_id = affs[0].get('comp_id')
                else:
                    for ent in (entry.get('nonpolymer_entities') or []):
                        try:
                            cid = ent['nonpolymer_comp']['chem_comp']['id'].upper()
                            if cid not in IGNORE: comp_id = cid; break
                        except: pass
            if not comp_id:
                self.error.emit('No ligand found for this PDB ID (apo or solvent only).'); return

            # 2. Download PDB
            pdb_resp = requests.get(f'https://files.rcsb.org/download/{self.pdb_id}.pdb', timeout=30)
            if not pdb_resp.ok:
                self.error.emit('Failed to download PDB from RCSB.'); return

            coords = []
            for line in pdb_resp.text.splitlines():
                if line.startswith('HETATM') and line[17:20].strip() == comp_id:
                    atom = line[12:16].strip()
                    if not atom.startswith('H'):
                        try:
                            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        except: pass

            if not coords:
                self.error.emit(f'Could not extract coordinates for ligand {comp_id}.'); return

            coords = np.array(coords)
            cx, cy, cz = coords.mean(axis=0)
            pad = 10.0
            sx = float(coords[:,0].max() - coords[:,0].min()) + pad
            sy = float(coords[:,1].max() - coords[:,1].min()) + pad
            sz = float(coords[:,2].max() - coords[:,2].min()) + pad
            self.finished.emit(cx, cy, cz, sx, sy, sz, comp_id)
        except Exception:
            self.error.emit(traceback.format_exc())


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN WINDOW
# ══════════════════════════════════════════════════════════════════════════════
class EasyDockingApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EASYDocking — Molecular Docking Platform")
        self.resize(1380, 900)
        self.setAcceptDrops(True)

        # State
        self._protein_data    = None   # str
        self._protein_path    = None   # temp file path
        self._protein_fmt     = 'pdbqt'
        self._ligands         = []     # list of (name, mol)
        self._results         = []     # list of dicts
        self._raw_pdb_bytes   = None
        self._ref_pdbqt       = None
        self._workers         = []
        self._current_result_row = 0

        central = QWidget(); self.setCentralWidget(central)
        root = QVBoxLayout(central); root.setContentsMargins(0,0,0,0); root.setSpacing(0)
        root.addWidget(self._make_header())

        self.tabs = QTabWidget()
        root.addWidget(self.tabs)
        self._build_protein_tab()
        self._build_ligand_tab()
        self._build_grid_tab()
        self._build_docking_tab()
        self._build_results_tab()
        self._build_interactions_tab()

        self.status = QStatusBar(); self.setStatusBar(self.status)
        self.status.showMessage("Ready — load a protein and ligands to start docking")

    # ── Header ─────────────────────────────────────────────────────────────────
    def _make_header(self):
        h = QFrame(); h.setFixedHeight(60)
        h.setStyleSheet(f"""QFrame{{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,
            stop:0 {C['bg']},stop:0.5 {C['surface']},stop:1 {C['bg']});
            border-bottom:1px solid {C['border2']};}}""")
        row = QHBoxLayout(h); row.setContentsMargins(20,0,20,0); row.setSpacing(12)

        icon = QLabel("🧬"); icon.setStyleSheet("font-size:26px;"); row.addWidget(icon)

        vlo = QVBoxLayout(); vlo.setContentsMargins(0,0,0,0); vlo.setSpacing(1)
        t = QLabel("EASYDocking")
        t.setFont(QFont("Inter", 17, QFont.Weight.Bold))
        t.setStyleSheet(f"color:{C['accent']};letter-spacing:0.5px;")
        sub = QLabel("Molecular Docking Platform  ·  AutoDock Vina")
        sub.setStyleSheet(f"color:{C['text2']};font-size:10px;")
        vlo.addWidget(t); vlo.addWidget(sub); row.addLayout(vlo)
        row.addStretch()

        for label, color in [("● Vina Engine", C['green']), ("● RDKit", C['green'])]:
            lbl = QLabel(label); lbl.setStyleSheet(f"color:{color};font-size:10px;")
            row.addWidget(lbl)

        badge = QLabel("v1.0  Desktop")
        badge.setStyleSheet(f"""color:{C['accent']};background:{C['surface2']};
            border:1px solid {C['border2']};border-radius:10px;
            padding:2px 10px;font-size:10px;font-weight:600;""")
        row.addWidget(badge)
        return h

    # ── Tab helpers ────────────────────────────────────────────────────────────
    def _section(self, title: str) -> QGroupBox:
        g = QGroupBox(title); return g

    def _primary_btn(self, label: str, color: str = None) -> QPushButton:
        c = color or C['accent']
        b = QPushButton(label)
        b.setStyleSheet(f"""QPushButton{{background:{c};color:#f0f4f9;font-weight:600;
            border:none;border-radius:6px;padding:7px 18px;}}
            QPushButton:hover{{background:{C['accent_h']};}}
            QPushButton:disabled{{background:{C['surface3']};color:{C['text3']};}}""")
        return b

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 1: PROTEIN
    # ══════════════════════════════════════════════════════════════════════════
    def _build_protein_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)

        # Input method
        inp_grp = self._section("Protein Input")
        inp_lo  = QVBoxLayout(inp_grp)

        mode_lo = QHBoxLayout()
        self.radio_prot_file = QRadioButton("Upload File  (.pdb / .pdbqt)")
        self.radio_prot_pdb  = QRadioButton("Fetch from RCSB  (PDB ID)")
        self.radio_prot_file.setChecked(True)
        bg = QButtonGroup(self)
        bg.addButton(self.radio_prot_file); bg.addButton(self.radio_prot_pdb)
        self.radio_prot_file.toggled.connect(self._toggle_prot_input)
        mode_lo.addWidget(self.radio_prot_file); mode_lo.addWidget(self.radio_prot_pdb)
        mode_lo.addStretch(); inp_lo.addLayout(mode_lo)

        # File upload panel
        self.prot_file_panel = QWidget()
        pfl = QHBoxLayout(self.prot_file_panel); pfl.setContentsMargins(0,0,0,0)
        self.btn_prot_browse = QPushButton("Browse Protein File…")
        self.btn_prot_browse.clicked.connect(self._browse_protein)
        self.lbl_prot_file = QLabel("No file loaded."); self.lbl_prot_file.setStyleSheet(f"color:{C['text2']};font-size:11px;")
        pfl.addWidget(self.btn_prot_browse); pfl.addWidget(self.lbl_prot_file); pfl.addStretch()
        inp_lo.addWidget(self.prot_file_panel)

        # PDB ID panel
        self.prot_pdb_panel = QWidget(); self.prot_pdb_panel.hide()
        pdbl = QHBoxLayout(self.prot_pdb_panel); pdbl.setContentsMargins(0,0,0,0)
        pdbl.addWidget(QLabel("PDB ID:"))
        self.inp_pdb_id = QLineEdit(); self.inp_pdb_id.setPlaceholderText("e.g. 1STP"); self.inp_pdb_id.setMaximumWidth(100)
        self.btn_fetch_pdb = self._primary_btn("Fetch from RCSB")
        self.btn_fetch_pdb.clicked.connect(self._fetch_pdb)
        pdbl.addWidget(self.inp_pdb_id); pdbl.addWidget(self.btn_fetch_pdb); pdbl.addStretch()
        inp_lo.addWidget(self.prot_pdb_panel)
        lo.addWidget(inp_grp)

        # Cleaning options
        clean_grp = self._section("Preprocessing")
        cl = QHBoxLayout(clean_grp)
        self.chk_rm_water  = QCheckBox("Remove Waters (HOH/WAT/SOL)"); self.chk_rm_water.setChecked(True)
        self.chk_rm_hetatm = QCheckBox("Remove HETATM (ligands, ions)"); self.chk_rm_hetatm.setChecked(True)
        cl.addWidget(self.chk_rm_water); cl.addWidget(self.chk_rm_hetatm); cl.addStretch()
        lo.addWidget(clean_grp)

        # Load button
        btn_lo = QHBoxLayout()
        self.btn_load_protein = self._primary_btn("⬆  Load & Prepare Protein")
        self.btn_load_protein.clicked.connect(self._load_protein)
        btn_lo.addWidget(self.btn_load_protein); btn_lo.addStretch()
        lo.addLayout(btn_lo)

        self.prot_console = Console()
        lo.addWidget(self.prot_console)

        # 3D viewer
        viewer_grp = self._section("Protein 3D Preview")
        vlo = QVBoxLayout(viewer_grp)
        self.prot_viewer = MolViewer()
        vlo.addWidget(self.prot_viewer)
        lo.addWidget(viewer_grp, 1)

        self.tabs.addTab(tab, "  🧫  Protein  ")

    def _toggle_prot_input(self):
        f = self.radio_prot_file.isChecked()
        self.prot_file_panel.setVisible(f)
        self.prot_pdb_panel.setVisible(not f)

    def _browse_protein(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Protein", "",
                                               "Protein Files (*.pdb *.pdbqt);;All Files (*)")
        if not path: return
        self.lbl_prot_file.setText(os.path.basename(path))
        with open(path, 'rb') as fh: self._raw_pdb_bytes = fh.read()
        self._pending_protein_path = path
        self._pending_protein_name = os.path.basename(path)

    def _fetch_pdb(self):
        pdb_id = self.inp_pdb_id.text().strip().upper()
        if len(pdb_id) != 4:
            QMessageBox.warning(self, "PDB ID", "PDB ID must be exactly 4 characters."); return
        self.prot_console.info(f"Fetching {pdb_id} from RCSB…")
        try:
            resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30)
            if not resp.ok:
                self.prot_console.err(f"Failed to fetch {pdb_id}: HTTP {resp.status_code}"); return
            self._raw_pdb_bytes = resp.content
            self._pending_protein_path = None
            self._pending_protein_name = f"{pdb_id}.pdb"
            self.prot_console.ok(f"{pdb_id} downloaded ({len(resp.content)//1024} KB)")
            self.lbl_prot_file.setText(f"{pdb_id}.pdb (from RCSB)")
        except Exception as e:
            self.prot_console.err(str(e))

    def _load_protein(self):
        if not self._raw_pdb_bytes:
            QMessageBox.warning(self, "No Protein", "Please select or fetch a protein file first."); return
        try:
            raw = self._raw_pdb_bytes.decode('utf-8', errors='replace')
            name = getattr(self, '_pending_protein_name', 'protein.pdb')

            # Clean
            lines = []
            for line in raw.splitlines():
                if self.chk_rm_hetatm.isChecked() and line.startswith('HETATM'): continue
                if self.chk_rm_water.isChecked() and len(line) >= 20:
                    if line[:6].strip() in ('ATOM','HETATM') and line[17:20].strip() in ('HOH','WAT','SOL'): continue
                lines.append(line)
            cleaned = '\n'.join(lines)

            if name.lower().endswith('.pdb'):
                self.prot_console.info("Converting PDB → PDBQT via OpenBabel…")
                pdbqt = self._pdb_to_pdbqt(cleaned)
                if not pdbqt:
                    self.prot_console.err("Conversion failed — check OpenBabel installation."); return
                self._protein_data = pdbqt
                self._protein_fmt  = 'pdbqt'
            else:
                self._protein_data = cleaned
                self._protein_fmt  = 'pdbqt'

            # Write temp file for Vina
            if self._protein_path and os.path.exists(self._protein_path):
                try: os.unlink(self._protein_path)
                except: pass
            tf = tempfile.NamedTemporaryFile(delete=False, suffix=f".{self._protein_fmt}")
            tf.write(self._protein_data.encode()); tf.close()
            self._protein_path = tf.name

            self.prot_console.ok(f"Protein ready: {name}  ({len(self._protein_data)//1024} KB)")
            self.prot_viewer.show_structure(self._protein_data, self._protein_fmt)
            self.status.showMessage(f"Protein loaded: {name}")

            # Pre-fill grid if PDB ID was used
            if hasattr(self, 'inp_pdb_id') and len(self.inp_pdb_id.text().strip()) == 4:
                self.inp_grid_pdb.setText(self.inp_pdb_id.text().strip().upper())

        except Exception as e:
            self.prot_console.err(traceback.format_exc())

    def _pdb_to_pdbqt(self, pdb_content: str) -> str | None:
        try:
            tin = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', mode='w')
            tin.write(pdb_content); tin.close()
            tout = tin.name + '.qt'
            r = subprocess.run([OBABEL,'-ipdb',tin.name,'-opdbqt','-O',tout,
                                 '-xr','-h','--partialcharge','gasteiger'],
                                capture_output=True, text=True)
            if r.returncode != 0:
                self.prot_console.warn(f"obabel stderr: {r.stderr[:200]}"); return None
            if os.path.exists(tout):
                with open(tout) as fh: content = fh.read()
                os.unlink(tout); os.unlink(tin.name)
                return content
            os.unlink(tin.name); return None
        except Exception as e:
            self.prot_console.err(str(e)); return None

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 2: LIGANDS
    # ══════════════════════════════════════════════════════════════════════════
    def _build_ligand_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)

        # Input mode tabs
        lig_tabs = QTabWidget()

        # ── SMILES ─────────────────────────────────────────────────────────
        smi_w = QWidget(); smi_lo = QVBoxLayout(smi_w); smi_lo.setContentsMargins(10,10,10,10)
        smi_lo.addWidget(QLabel("Paste SMILES — one per line. Format: SMILES  Name  (Name optional)",
                                styleSheet=f"color:{C['text2']};font-size:11px;"))
        self.txt_smiles = QTextEdit(); self.txt_smiles.setPlaceholderText(
            "CCO  Ethanol\nC1=CC=CC=C1  Benzene\n...")
        smi_lo.addWidget(self.txt_smiles)
        lig_tabs.addTab(smi_w, "📝  SMILES")

        # ── SDF ────────────────────────────────────────────────────────────
        sdf_w = QWidget(); sdf_lo = QVBoxLayout(sdf_w); sdf_lo.setContentsMargins(10,10,10,10)
        sdf_btn_lo = QHBoxLayout()
        self.btn_load_sdf = QPushButton("Browse SDF File…"); self.btn_load_sdf.clicked.connect(self._load_sdf)
        self.lbl_sdf = QLabel("No SDF loaded."); self.lbl_sdf.setStyleSheet(f"color:{C['text2']};font-size:11px;")
        sdf_btn_lo.addWidget(self.btn_load_sdf); sdf_btn_lo.addWidget(self.lbl_sdf); sdf_btn_lo.addStretch()
        sdf_lo.addLayout(sdf_btn_lo)
        lig_tabs.addTab(sdf_w, "📂  SDF File")

        # ── Redocking ──────────────────────────────────────────────────────
        red_w = QWidget(); red_lo = QVBoxLayout(red_w); red_lo.setContentsMargins(10,10,10,10)
        red_lo.addWidget(QLabel("Extract co-crystallized ligand from the loaded protein PDB.",
                                styleSheet=f"color:{C['text2']};font-size:11px;"))
        rl = QHBoxLayout()
        rl.addWidget(QLabel("Residue code (3 letters):"))
        self.inp_resname = QLineEdit(); self.inp_resname.setMaximumWidth(80)
        self.inp_resname.setPlaceholderText("e.g. LIG")
        self.btn_extract = QPushButton("Extract Ligand"); self.btn_extract.clicked.connect(self._extract_ligand)
        self.lbl_extract = QLabel(""); self.lbl_extract.setStyleSheet(f"color:{C['text2']};font-size:11px;")
        rl.addWidget(self.inp_resname); rl.addWidget(self.btn_extract); rl.addWidget(self.lbl_extract); rl.addStretch()
        red_lo.addLayout(rl)
        lig_tabs.addTab(red_w, "♻️  Redocking")

        lo.addWidget(lig_tabs)

        # Advanced options
        adv_grp = self._section("Advanced Ligand Preparation")
        adv_lo  = QHBoxLayout(adv_grp)
        self.chk_minimize = QCheckBox("Minimize geometry (MMFF94) — recommended for redocking")
        adv_lo.addWidget(self.chk_minimize); adv_lo.addStretch()
        lo.addWidget(adv_grp)

        # Load button
        btn_lo = QHBoxLayout()
        self.btn_load_ligands = self._primary_btn("⬆  Load Ligands", C['purple'])
        self.btn_load_ligands.clicked.connect(self._load_ligands)
        btn_lo.addWidget(self.btn_load_ligands); btn_lo.addStretch()
        self.lbl_lig_count = QLabel("0 ligands loaded")
        self.lbl_lig_count.setStyleSheet(f"color:{C['text2']};font-size:11px;")
        btn_lo.addWidget(self.lbl_lig_count)
        lo.addLayout(btn_lo)

        self.lig_console = Console()
        lo.addWidget(self.lig_console)

        # Ligand table
        tbl_grp = self._section("Loaded Ligands")
        tbl_lo  = QVBoxLayout(tbl_grp)
        self.tbl_ligands = QTableWidget(0, 3)
        self.tbl_ligands.setHorizontalHeaderLabels(['Name','SMILES / Source','Heavy Atoms'])
        self.tbl_ligands.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_ligands.setMaximumHeight(200)
        tbl_lo.addWidget(self.tbl_ligands)
        lo.addWidget(tbl_grp)
        self.tabs.addTab(tab, "  💊  Ligands  ")

    def _load_sdf(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load SDF", "", "SDF Files (*.sdf);;All Files (*)")
        if not path: return
        self._pending_sdf_path = path
        self.lbl_sdf.setText(os.path.basename(path))

    def _extract_ligand(self):
        if not self._raw_pdb_bytes:
            QMessageBox.warning(self, "No Protein", "Load a protein first."); return
        resname = self.inp_resname.text().strip().upper()
        if len(resname) != 3:
            QMessageBox.warning(self, "Residue Code", "Enter a 3-letter residue code."); return
        raw = self._raw_pdb_bytes.decode('utf-8', errors='replace')
        lines = [l for l in raw.splitlines()
                 if l[:6] in ('HETATM','ATOM  ') and l[17:20].strip() == resname]
        if not lines:
            self.lbl_extract.setText(f"⚠ Residue {resname} not found"); return
        self._pending_redock_lines = '\n'.join(lines)
        self._pending_redock_name  = resname
        self.lbl_extract.setText(f"✔ {len(lines)} atoms found for {resname}")

    def _load_ligands(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        self._ligands = []

        # SMILES
        txt = self.txt_smiles.toPlainText().strip()
        if txt:
            for i, line in enumerate(txt.splitlines()):
                parts = line.strip().split()
                if not parts: continue
                smi  = parts[0]
                name = parts[1] if len(parts) > 1 else f"Ligand_{i+1}"
                mol  = Chem.MolFromSmiles(smi)
                if mol:
                    mol = Chem.AddHs(mol)
                    self._ligands.append((name, mol))
                else:
                    self.lig_console.warn(f"Invalid SMILES: {smi}")

        # SDF
        if hasattr(self, '_pending_sdf_path') and self._pending_sdf_path:
            suppl = Chem.SDMolSupplier(self._pending_sdf_path)
            for i, mol in enumerate(suppl):
                if mol:
                    name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"SDF_{i+1}"
                    mol  = Chem.AddHs(mol)
                    self._ligands.append((name, mol))

        # Redocking
        if hasattr(self, '_pending_redock_lines') and self._pending_redock_lines:
            resname = self._pending_redock_name
            try:
                with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as t:
                    t.write(self._pending_redock_lines); tpath = t.name
                tout = tpath.replace('.pdb', '.sdf')
                subprocess.run([OBABEL,'-ipdb',tpath,'-osdf','-O',tout,'-h'],
                               capture_output=True)
                if os.path.exists(tout):
                    suppl = Chem.SDMolSupplier(tout, removeHs=False)
                    mols  = [m for m in suppl if m]
                    if mols:
                        mol = Chem.AddHs(mols[0])
                        self._ligands.append((f"{resname}_Redocking", mol))
                        # Save ref PDBQT for visualization
                        tpdbqt = tpath.replace('.pdb', '.pdbqt')
                        subprocess.run([OBABEL,'-ipdb',tpath,'-opdbqt','-O',tpdbqt,
                                        '-h','--partialcharge','gasteiger'], capture_output=True)
                        if os.path.exists(tpdbqt):
                            with open(tpdbqt) as fh: self._ref_pdbqt = fh.read()
                            os.unlink(tpdbqt)
                    os.unlink(tout)
                os.unlink(tpath)
            except Exception as e:
                self.lig_console.err(str(e))

        # Populate table
        self.tbl_ligands.setRowCount(len(self._ligands))
        for i, (name, mol) in enumerate(self._ligands):
            from rdkit import Chem
            smi = Chem.MolToSmiles(mol)
            for j, val in enumerate([name, smi[:60]+'…' if len(smi)>60 else smi,
                                      str(mol.GetNumHeavyAtoms())]):
                self.tbl_ligands.setItem(i, j, QTableWidgetItem(val))
        n = len(self._ligands)
        self.lbl_lig_count.setText(f"{n} ligand{'s' if n!=1 else ''} loaded")
        self.lig_console.ok(f"{n} ligands ready")
        if n > 0: self.status.showMessage(f"{n} ligands loaded")

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 3: GRID
    # ══════════════════════════════════════════════════════════════════════════
    def _build_grid_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left: controls
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0,0,0,0); ll.setSpacing(10)

        # Auto from PDB
        auto_grp = self._section("Auto-Calculate from PDB ID")
        auto_lo  = QHBoxLayout(auto_grp)
        self.inp_grid_pdb = QLineEdit(); self.inp_grid_pdb.setPlaceholderText("e.g. 1STP"); self.inp_grid_pdb.setMaximumWidth(100)
        self.btn_calc_grid = self._primary_btn("Calculate Grid")
        self.btn_calc_grid.clicked.connect(self._calc_grid_auto)
        self.lbl_grid_info = QLabel(""); self.lbl_grid_info.setStyleSheet(f"color:{C['green']};font-size:10px;")
        auto_lo.addWidget(QLabel("PDB ID:")); auto_lo.addWidget(self.inp_grid_pdb)
        auto_lo.addWidget(self.btn_calc_grid); auto_lo.addWidget(self.lbl_grid_info); auto_lo.addStretch()
        ll.addWidget(auto_grp)

        # Manual coordinates
        coord_grp = self._section("Grid Box Coordinates")
        cg = QVBoxLayout(coord_grp); cg.setSpacing(6)

        def _spin(val=0.0, rng=(-999,999)):
            s = QDoubleSpinBox(); s.setRange(*rng); s.setValue(val); s.setSingleStep(0.5)
            s.setDecimals(2); s.setMinimumWidth(100); return s

        def _row(label, sp):
            h = QHBoxLayout()
            l = QLabel(label); l.setFixedWidth(90); l.setStyleSheet(f"color:{C['text2']};")
            h.addWidget(l); h.addWidget(sp); h.addStretch(); return h

        self.spin_cx = _spin(); self.spin_cy = _spin(); self.spin_cz = _spin()
        self.spin_sx = _spin(20.0, (1,200)); self.spin_sy = _spin(20.0, (1,200)); self.spin_sz = _spin(20.0, (1,200))

        for lbl, sp in [("Center X:", self.spin_cx), ("Center Y:", self.spin_cy), ("Center Z:", self.spin_cz),
                         ("Size X (Å):", self.spin_sx), ("Size Y (Å):", self.spin_sy), ("Size Z (Å):", self.spin_sz)]:
            cg.addLayout(_row(lbl, sp))

        for sp in [self.spin_cx, self.spin_cy, self.spin_cz,
                   self.spin_sx, self.spin_sy, self.spin_sz]:
            sp.valueChanged.connect(self._update_grid_preview)

        ll.addWidget(coord_grp)

        # Preview button
        btn_lo = QHBoxLayout()
        btn_preview = self._primary_btn("🔄  Update 3D Grid Preview")
        btn_preview.clicked.connect(self._update_grid_preview)
        btn_lo.addWidget(btn_preview); btn_lo.addStretch()
        ll.addLayout(btn_lo)

        self.grid_console = Console(); self.grid_console.setMaximumHeight(100)
        ll.addWidget(self.grid_console)
        ll.addStretch()
        splitter.addWidget(left)

        # Right: 3D viewer
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0,0,0,0)
        lbl = QLabel("3D Preview — Protein + Grid Box", styleSheet=f"color:{C['text2']};font-size:10px;")
        rl.addWidget(lbl)
        self.grid_viewer = MolViewer(); rl.addWidget(self.grid_viewer)
        splitter.addWidget(right)
        splitter.setSizes([380, 900])

        lo.addWidget(splitter, 1)
        self.tabs.addTab(tab, "  📦  Grid Box  ")

    def _calc_grid_auto(self):
        pdb_id = self.inp_grid_pdb.text().strip().upper()
        if len(pdb_id) != 4:
            QMessageBox.warning(self, "PDB ID", "Enter a 4-character PDB ID."); return
        self.grid_console.info(f"Fetching grid parameters for {pdb_id}…")
        self.btn_calc_grid.setEnabled(False)
        w = GridWorker(pdb_id)
        w.finished.connect(self._on_grid_done)
        w.error.connect(lambda e: (self.grid_console.err(e), self.btn_calc_grid.setEnabled(True)))
        self._workers.append(w); w.start()

    def _on_grid_done(self, cx, cy, cz, sx, sy, sz, comp_id):
        self.btn_calc_grid.setEnabled(True)
        for sp, v in [(self.spin_cx,cx),(self.spin_cy,cy),(self.spin_cz,cz),
                      (self.spin_sx,sx),(self.spin_sy,sy),(self.spin_sz,sz)]:
            sp.blockSignals(True); sp.setValue(v); sp.blockSignals(False)
        self.lbl_grid_info.setText(f"✔ Grid calculated from ligand {comp_id}")
        self.grid_console.ok(f"Grid from {comp_id}: center ({cx:.1f}, {cy:.1f}, {cz:.1f}) size ({sx:.1f}, {sy:.1f}, {sz:.1f})")
        self._update_grid_preview()

    def _update_grid_preview(self):
        if not self._protein_data: return
        cx=self.spin_cx.value(); cy=self.spin_cy.value(); cz=self.spin_cz.value()
        sx=self.spin_sx.value(); sy=self.spin_sy.value(); sz=self.spin_sz.value()
        self.grid_viewer.show_structure(self._protein_data, self._protein_fmt,
                                        center=(cx,cy,cz), box_size=(sx,sy,sz))

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 4: DOCKING
    # ══════════════════════════════════════════════════════════════════════════
    def _build_docking_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)

        # Parameters
        param_grp = self._section("Docking Parameters")
        pg = QHBoxLayout(param_grp)

        def _lbl_spin(label, sp):
            h = QHBoxLayout()
            l = QLabel(label); l.setStyleSheet(f"color:{C['text2']};")
            h.addWidget(l); h.addWidget(sp); h.addStretch()
            return h

        self.spin_exh = QSpinBox(); self.spin_exh.setRange(1,32); self.spin_exh.setValue(8); self.spin_exh.setMinimumWidth(80)
        self.spin_poses = QSpinBox(); self.spin_poses.setRange(1,20); self.spin_poses.setValue(9); self.spin_poses.setMinimumWidth(80)
        pg.addLayout(_lbl_spin("Exhaustiveness:", self.spin_exh))
        pg.addLayout(_lbl_spin("Num Poses:", self.spin_poses))
        pg.addStretch()
        lo.addWidget(param_grp)

        # Summary
        sum_grp = self._section("Ready to Dock?")
        sg = QVBoxLayout(sum_grp)
        self.lbl_dock_summary = QLabel("—")
        self.lbl_dock_summary.setStyleSheet(f"color:{C['text2']};font-size:11px;")
        sg.addWidget(self.lbl_dock_summary)
        lo.addWidget(sum_grp)

        # Run button
        btn_lo = QHBoxLayout()
        self.btn_run_dock = self._primary_btn("▶   Run Docking", C['green'])
        self.btn_run_dock.setFixedHeight(46)
        self.btn_run_dock.clicked.connect(self._run_docking)
        btn_lo.addWidget(self.btn_run_dock, 1)
        lo.addLayout(btn_lo)

        # Progress
        self.dock_progress = QProgressBar(); self.dock_progress.setValue(0)
        lo.addWidget(self.dock_progress)

        self.dock_console = Console(); self.dock_console.setMaximumHeight(220)
        lo.addWidget(self.dock_console)
        lo.addStretch()
        self.tabs.addTab(tab, "  ⚗️  Docking  ")

    def _run_docking(self):
        # Validate
        msgs = []
        if not self._protein_path or not os.path.exists(self._protein_path):
            msgs.append("No protein loaded (Tab 1)")
        if not self._ligands:
            msgs.append("No ligands loaded (Tab 2)")
        if msgs:
            QMessageBox.warning(self, "Not Ready", '\n'.join(msgs)); return

        cx=self.spin_cx.value(); cy=self.spin_cy.value(); cz=self.spin_cz.value()
        sx=self.spin_sx.value(); sy=self.spin_sy.value(); sz=self.spin_sz.value()

        params = {'cx':cx,'cy':cy,'cz':cz,'sx':sx,'sy':sy,'sz':sz,
                  'exhaustiveness': self.spin_exh.value(),
                  'n_poses':        self.spin_poses.value(),
                  'minimize':       self.chk_minimize.isChecked() if hasattr(self,'chk_minimize') else False}

        self._results = []
        self.btn_run_dock.setEnabled(False)
        self.dock_progress.setValue(0)
        self.dock_console.info(f"Starting docking of {len(self._ligands)} ligands…")
        self.dock_console.info(f"Grid: center ({cx:.1f}, {cy:.1f}, {cz:.1f})  size ({sx:.1f}, {sy:.1f}, {sz:.1f})")

        w = DockingWorker(self._protein_path, self._ligands, params)
        w.progress.connect(self._on_dock_progress)
        w.result.connect(self._on_dock_result)
        w.finished.connect(self._on_dock_finished)
        w.error.connect(lambda e: (self.dock_console.err(e), self.btn_run_dock.setEnabled(True)))
        self._workers.append(w); w.start()

    def _on_dock_progress(self, done, total, name):
        pct = int((done / max(total,1)) * 100)
        self.dock_progress.setValue(pct)
        self.dock_console.info(f"[{done+1}/{total}] Docking {name}…")

    def _on_dock_result(self, r: dict):
        self._results.append(r)
        if 'error' in r:
            self.dock_console.warn(f"{r['Ligand']}: {r['error']}")
        else:
            self.dock_console.ok(f"{r['Ligand']}  Affinity: {r['Affinity (kcal/mol)']} kcal/mol")

    def _on_dock_finished(self):
        self.dock_progress.setValue(100)
        self.btn_run_dock.setEnabled(True)
        self.dock_console.ok(f"Docking complete — {len(self._results)} results")
        self.status.showMessage(f"Docking done: {len(self._results)} ligands processed")
        self._populate_results()
        self._populate_interaction_combo()
        self.tabs.setCurrentIndex(4)

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 5: RESULTS
    # ══════════════════════════════════════════════════════════════════════════
    def _build_results_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)

        splitter = QSplitter(Qt.Orientation.Vertical)

        # Top: table + export
        top = QWidget(); tl = QVBoxLayout(top); tl.setContentsMargins(0,0,0,0)

        tbl_grp = self._section("Docking Results")
        tbl_vlo = QVBoxLayout(tbl_grp)
        self.tbl_results = QTableWidget(0, 5)
        self.tbl_results.setHorizontalHeaderLabels(
            ['Ligand','Affinity (kcal/mol)','RMSD (Å)','Heavy Atoms','Ligand Efficiency'])
        self.tbl_results.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_results.setAlternatingRowColors(True)
        self.tbl_results.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tbl_results.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.tbl_results.itemSelectionChanged.connect(self._on_result_selected)
        tbl_vlo.addWidget(self.tbl_results)

        # Export buttons
        exp_lo = QHBoxLayout()
        self.btn_exp_csv = self._primary_btn("📥 Export CSV")
        self.btn_exp_zip = self._primary_btn("📦 Export ZIP (All Files)", C['purple'])
        self.btn_exp_csv.clicked.connect(self._export_csv)
        self.btn_exp_zip.clicked.connect(self._export_zip)
        exp_lo.addWidget(self.btn_exp_csv); exp_lo.addWidget(self.btn_exp_zip); exp_lo.addStretch()
        tbl_vlo.addLayout(exp_lo)
        tl.addWidget(tbl_grp)
        splitter.addWidget(top)

        # Bottom: 3D viewer + affinity chart
        bottom = QWidget(); bl = QHBoxLayout(bottom); bl.setContentsMargins(0,0,0,0); bl.setSpacing(8)

        viz_grp = self._section("3D Pose Viewer")
        vl = QVBoxLayout(viz_grp)
        self.result_viewer = MolViewer(); self.result_viewer.setMinimumHeight(340)
        self.lbl_viz_info = QLabel("Select a ligand from the table above to visualize its best pose.",
                                   styleSheet=f"color:{C['text2']};font-size:10px;")
        vl.addWidget(self.result_viewer); vl.addWidget(self.lbl_viz_info)
        bl.addWidget(viz_grp, 2)

        chart_grp = self._section("Affinity Chart")
        cl = QVBoxLayout(chart_grp)
        self.result_fig = Figure(facecolor=C['bg'])
        self.result_canvas = FigureCanvas(self.result_fig)
        cl.addWidget(self.result_canvas)
        bl.addWidget(chart_grp, 1)

        splitter.addWidget(bottom)
        splitter.setSizes([280, 420])
        lo.addWidget(splitter, 1)
        self.tabs.addTab(tab, "  📊  Results  ")

    def _populate_results(self):
        good = [r for r in self._results if 'error' not in r]
        good.sort(key=lambda x: x.get('Affinity (kcal/mol)', 0))

        cols = ['Ligand','Affinity (kcal/mol)','RMSD (Å)','Heavy Atoms','Ligand Efficiency']
        self.tbl_results.setRowCount(len(good))
        for i, r in enumerate(good):
            for j, col in enumerate(cols):
                val = str(r.get(col, '—'))
                item = QTableWidgetItem(val)
                if col == 'Affinity (kcal/mol)':
                    try:
                        v = float(val)
                        item.setForeground(QColor(C['green'] if v < -8 else
                                                   C['amber'] if v < -6 else C['red']))
                    except: pass
                self.tbl_results.setItem(i, j, item)

        # Select best
        if good: self.tbl_results.selectRow(0)
        self._draw_affinity_chart(good)

    def _on_result_selected(self):
        rows = self.tbl_results.selectedItems()
        if not rows or not self._protein_data: return
        row = rows[0].row()
        lig_item = self.tbl_results.item(row, 0)
        if not lig_item: return
        lig_name = lig_item.text()

        r = next((x for x in self._results if x.get('Ligand') == lig_name), None)
        if not r or 'PDBQT' not in r: return

        cx=self.spin_cx.value(); cy=self.spin_cy.value(); cz=self.spin_cz.value()
        sx=self.spin_sx.value(); sy=self.spin_sy.value(); sz=self.spin_sz.value()

        ref = self._ref_pdbqt if '_Redocking' in lig_name else None
        self.result_viewer.show_structure(
            self._protein_data, self._protein_fmt,
            ligand_pdbqt=r['PDBQT'],
            center=(cx,cy,cz), box_size=(sx,sy,sz),
            ref_pdbqt=ref)
        aff = r.get('Affinity (kcal/mol)', '—')
        self.lbl_viz_info.setText(
            f"Ligand: {lig_name}  |  Affinity: {aff} kcal/mol  "
            + ("  |  Cyan = native pose" if ref else ""))

    def _draw_affinity_chart(self, results: list):
        self.result_fig.clear()
        ax = self.result_fig.add_subplot(111)
        ax.set_facecolor(C['surface2'])
        ax.tick_params(colors=C['text2'], labelsize=8)
        for sp in ax.spines.values(): sp.set_color(C['border'])
        ax.grid(True, color=C['border'], linewidth=0.4, alpha=0.6, axis='y')

        if not results: return
        names = [r['Ligand'][:20]+'…' if len(r['Ligand'])>20 else r['Ligand'] for r in results]
        vals  = [r.get('Affinity (kcal/mol)', 0) for r in results]
        colors_bar = [C['green'] if v < -8 else C['amber'] if v < -6 else C['red'] for v in vals]

        bars = ax.bar(range(len(names)), vals, color=colors_bar, alpha=0.8)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=45, ha='right', fontsize=7)
        ax.set_ylabel('Affinity (kcal/mol)', color=C['text2'], fontsize=8)
        ax.set_title('Binding Affinities', color=C['text'], fontsize=9)
        ax.axhline(-8, color=C['green'], ls='--', lw=0.8, alpha=0.6)
        ax.axhline(-6, color=C['amber'], ls='--', lw=0.8, alpha=0.6)
        self.result_fig.tight_layout()
        self.result_canvas.draw()

    def _export_csv(self):
        good = [r for r in self._results if 'error' not in r]
        if not good: QMessageBox.information(self,"No Results","Run docking first."); return
        path, _ = QFileDialog.getSaveFileName(self,"Save CSV","docking_results.csv","CSV (*.csv)")
        if not path: return
        cols = ['Ligand','Affinity (kcal/mol)','RMSD (Å)','Heavy Atoms','Ligand Efficiency']
        df = pd.DataFrame([{c: r.get(c,'') for c in cols} for r in good])
        df.to_csv(path, index=False)
        QMessageBox.information(self,"Exported",f"✔ Saved to:\n{path}")

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 6: INTERACTIONS
    # ══════════════════════════════════════════════════════════════════════════
    def _build_interactions_tab(self):
        tab = QWidget(); lo = QVBoxLayout(tab); lo.setContentsMargins(14,14,14,14); lo.setSpacing(10)

        # Top controls
        ctrl_lo = QHBoxLayout()
        ctrl_lo.addWidget(QLabel("Ligand:", styleSheet=f"color:{C['text2']};"))
        self.combo_ilig = QComboBox(); self.combo_ilig.setMinimumWidth(260)
        ctrl_lo.addWidget(self.combo_ilig)
        self.btn_run_int = self._primary_btn("🔬  Analyze Interactions", C['purple'])
        self.btn_run_int.clicked.connect(self._run_interactions)
        self.btn_run_int.setEnabled(False)
        ctrl_lo.addWidget(self.btn_run_int)
        ctrl_lo.addStretch()
        self.lbl_int_status = QLabel("Run docking first, then select a ligand.")
        self.lbl_int_status.setStyleSheet(f"color:{C['text2']};font-size:10px;")
        ctrl_lo.addWidget(self.lbl_int_status)
        lo.addLayout(ctrl_lo)

        self.int_pbar = QProgressBar(); self.int_pbar.setRange(0, 0); self.int_pbar.setFixedHeight(4)
        self.int_pbar.hide(); lo.addWidget(self.int_pbar)

        # KPI cards row
        kpi_frame = QFrame()
        kpi_frame.setStyleSheet(f"QFrame{{background:{C['surface']};border-radius:8px;}}")
        kpi_lo = QHBoxLayout(kpi_frame); kpi_lo.setContentsMargins(12,8,12,8); kpi_lo.setSpacing(0)
        self._int_kpis = {}
        kpi_defs = [
            ('Total Residues',  C['accent']),
            ('H-Bonds',         '#3d8ef0'),
            ('Hydrophobic',     C['amber']),
            ('π–π / π–Cat',    C['purple']),
            ('Ionic',           C['red']),
            ('VdW Contacts',    C['text2']),
        ]
        for i, (name, color) in enumerate(kpi_defs):
            vlo = QVBoxLayout(); vlo.setSpacing(2)
            val = QLabel('—'); val.setFont(QFont('Inter', 18, QFont.Weight.Bold))
            val.setStyleSheet(f'color:{color};'); val.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl = QLabel(name); lbl.setStyleSheet(f'color:{C["text2"]};font-size:9px;font-weight:700;letter-spacing:0.5px;')
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            vlo.addWidget(val); vlo.addWidget(lbl); kpi_lo.addLayout(vlo)
            self._int_kpis[name] = val
            if i < len(kpi_defs) - 1:
                sep = QFrame(); sep.setFrameShape(QFrame.Shape.VLine)
                sep.setStyleSheet(f'color:{C["border"]};margin:10px 16px;'); kpi_lo.addWidget(sep)
        lo.addWidget(kpi_frame)

        # Main split: diagram (left) + table (right)
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left: 2D diagram
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0,0,0,0)
        ll.addWidget(QLabel("2D Interaction Diagram",
                            styleSheet=f"color:{C['text2']};font-size:10px;font-weight:700;"))
        self.int_diagram = InteractionDiagram()
        ll.addWidget(self.int_diagram, 1)

        # Export diagram
        btn_exp_diag = QPushButton("💾 Save Diagram (PNG)")
        btn_exp_diag.setStyleSheet("padding:5px 12px;")
        btn_exp_diag.clicked.connect(self._export_diagram)
        ll.addWidget(btn_exp_diag)
        splitter.addWidget(left)

        # Right: interaction table + type chart
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0,0,0,0); rl.setSpacing(8)

        # Table
        tbl_lbl = QLabel("Interaction Detail Table",
                          styleSheet=f"color:{C['text2']};font-size:10px;font-weight:700;")
        rl.addWidget(tbl_lbl)
        self.tbl_int = QTableWidget(0, 3)
        self.tbl_int.setHorizontalHeaderLabels(['Residue', 'Interaction Type', '●'])
        self.tbl_int.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tbl_int.setAlternatingRowColors(True)
        self.tbl_int.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tbl_int.setMaximumHeight(260)
        rl.addWidget(self.tbl_int)

        # Bar chart: interactions per type
        chart_lbl = QLabel("Interaction Count by Type",
                            styleSheet=f"color:{C['text2']};font-size:10px;font-weight:700;")
        rl.addWidget(chart_lbl)
        self.int_fig = Figure(facecolor=C['bg'], figsize=(5, 3))
        self.int_canvas = FigureCanvas(self.int_fig); self.int_canvas.setMinimumHeight(200)
        rl.addWidget(self.int_canvas, 1)

        # Export CSV
        btn_exp_csv = QPushButton("📥 Export Interactions (CSV)")
        btn_exp_csv.setStyleSheet("padding:5px 12px;")
        btn_exp_csv.clicked.connect(self._export_interactions_csv)
        rl.addWidget(btn_exp_csv)

        splitter.addWidget(right)
        splitter.setSizes([680, 440])
        lo.addWidget(splitter, 1)

        self._last_interactions = {}
        self._last_int_ligand   = ''
        self.tabs.addTab(tab, "  🔬  Interactions  ")

    # ── Interaction methods ────────────────────────────────────────────────────

    def _populate_interaction_combo(self):
        good = [r for r in self._results if 'error' not in r and 'PDBQT' in r]
        self.combo_ilig.clear()
        for r in good:
            self.combo_ilig.addItem(r['Ligand'])
        if good:
            self.btn_run_int.setEnabled(True)
            # Auto-run on best result
            self.tabs.setCurrentIndex(5)
            self._run_interactions()

    def _run_interactions(self):
        if not self._protein_data:
            QMessageBox.warning(self, "No Protein", "Load a protein first."); return
        lig_name = self.combo_ilig.currentText()
        if not lig_name:
            QMessageBox.warning(self, "No Ligand", "No docking results available."); return

        r = next((x for x in self._results
                  if x.get('Ligand') == lig_name and 'PDBQT' in x), None)
        if not r:
            QMessageBox.warning(self, "No Result", f"No PDBQT found for {lig_name}."); return

        self.btn_run_int.setEnabled(False)
        self.int_pbar.show()
        self.lbl_int_status.setText(f"Analyzing {lig_name}…")
        self.lbl_int_status.setStyleSheet(f"color:{C['text2']};font-size:10px;")
        self.int_diagram._placeholder()

        # Get ligand SMILES
        lig_smi = ''
        for name, mol in self._ligands:
            if name == lig_name:
                try:
                    from rdkit import Chem
                    lig_smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
                except: pass
                break

        w = InteractionWorker(self._protein_data, r['PDBQT'], lig_name)
        w.progress.connect(self.lbl_int_status.setText)
        w.finished.connect(self._on_interactions_done)
        w.error.connect(lambda e: (
            self.lbl_int_status.setText(f"⚠ Error — check console"),
            self.lbl_int_status.setStyleSheet(f"color:{C['red']};font-size:10px;"),
            self.int_pbar.hide(),
            self.btn_run_int.setEnabled(True)
        ))
        self._workers.append(w); w.start()

    def _on_interactions_done(self, data: dict):
        self.int_pbar.hide()
        self.btn_run_int.setEnabled(True)

        interactions = data.get('interactions', {})
        lig_name     = data.get('ligand_name', '')
        lig_smi      = data.get('ligand_smi', '')
        self._last_interactions = interactions
        self._last_int_ligand   = lig_name

        if not interactions:
            self.lbl_int_status.setText("No interactions found — try a different pose or check grid.")
            self.lbl_int_status.setStyleSheet(f"color:{C['amber']};font-size:10px;")
            return

        # ── KPIs ──────────────────────────────────────────────────────────────
        all_types = [t for ts in interactions.values() for t in ts]
        hbonds    = sum(1 for t in all_types if 'HB' in t)
        hydro     = sum(1 for t in all_types if 'Hydrophobic' in t)
        pi        = sum(1 for t in all_types if 'Pi' in t or 'Cation' in t)
        ionic     = sum(1 for t in all_types if t in ('Cationic','Anionic'))
        vdw       = sum(1 for t in all_types if t == 'VdWContact')

        self._int_kpis['Total Residues'].setText(str(len(interactions)))
        self._int_kpis['H-Bonds'].setText(str(hbonds))
        self._int_kpis['Hydrophobic'].setText(str(hydro))
        self._int_kpis['π–π / π–Cat'].setText(str(pi))
        self._int_kpis['Ionic'].setText(str(ionic))
        self._int_kpis['VdW Contacts'].setText(str(vdw))

        # ── Table ──────────────────────────────────────────────────────────────
        rows = [(res, itype)
                for res, itypes in interactions.items()
                for itype in itypes]
        self.tbl_int.setRowCount(len(rows))
        for i, (res, itype) in enumerate(rows):
            style = INTERACTION_STYLES.get(itype, INTERACTION_STYLES['VdWContact'])
            self.tbl_int.setItem(i, 0, QTableWidgetItem(str(res)))
            itype_item = QTableWidgetItem(style['label'])
            itype_item.setForeground(QColor(style['color']))
            self.tbl_int.setItem(i, 1, itype_item)
            dot_item = QTableWidgetItem('●')
            dot_item.setForeground(QColor(style['color']))
            self.tbl_int.setItem(i, 2, dot_item)

        # ── Bar chart ──────────────────────────────────────────────────────────
        type_counts = {}
        for itype in all_types:
            st = INTERACTION_STYLES.get(itype, INTERACTION_STYLES['VdWContact'])
            lbl = st['label']
            type_counts[lbl] = type_counts.get(lbl, 0) + 1

        self.int_fig.clear()
        ax = self.int_fig.add_subplot(111)
        ax.set_facecolor(C['surface2'])
        ax.tick_params(colors=C['text2'], labelsize=8)
        for sp in ax.spines.values(): sp.set_color(C['border'])
        ax.grid(True, color=C['border'], linewidth=0.4, alpha=0.5, axis='x')

        sorted_types = sorted(type_counts.items(), key=lambda x: x[1], reverse=True)
        labels = [x[0] for x in sorted_types]
        values = [x[1] for x in sorted_types]
        # Match colors
        bar_colors = []
        for lbl in labels:
            c_val = next((st['color'] for st in INTERACTION_STYLES.values()
                         if st['label'] == lbl), C['text2'])
            bar_colors.append(c_val)

        bars = ax.barh(range(len(labels)), values, color=bar_colors, alpha=0.82)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Count', color=C['text2'], fontsize=8)
        ax.set_title('Interactions by Type', color=C['text'], fontsize=9)
        for bar, val in zip(bars, values):
            ax.text(val + 0.05, bar.get_y() + bar.get_height()/2,
                    str(val), va='center', color=C['text2'], fontsize=8)
        self.int_fig.tight_layout()
        self.int_canvas.draw()

        # ── 2D diagram ─────────────────────────────────────────────────────────
        self.int_diagram.render(interactions, lig_name, lig_smi)

        n = len(interactions)
        self.lbl_int_status.setText(f"✔ {n} interacting residues  |  {len(all_types)} total interactions")
        self.lbl_int_status.setStyleSheet(f"color:{C['green']};font-size:10px;")
        self.status.showMessage(f"Interactions: {n} residues — {lig_name}")

    def _export_diagram(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Diagram", f"interactions_{self._last_int_ligand}.png",
            "PNG Image (*.png);;PDF (*.pdf)")
        if path:
            self.int_diagram.fig.savefig(path, dpi=180, bbox_inches='tight',
                                          facecolor=C['bg'])
            QMessageBox.information(self, "Saved", f"Diagram saved to:\n{path}")

    def _export_interactions_csv(self):
        if not self._last_interactions:
            QMessageBox.information(self, "No Data", "Run interaction analysis first."); return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV", f"interactions_{self._last_int_ligand}.csv",
            "CSV Files (*.csv)")
        if not path: return
        rows = []
        for res, itypes in self._last_interactions.items():
            for itype in itypes:
                rows.append({'Ligand': self._last_int_ligand,
                             'Residue': res, 'Interaction': itype,
                             'Category': INTERACTION_STYLES.get(itype, {}).get('label', itype)})
        pd.DataFrame(rows).to_csv(path, index=False)
        QMessageBox.information(self, "Exported", f"✔ Saved to:\n{path}")

    def _export_zip(self):
        good = [r for r in self._results if 'error' not in r]
        if not good: QMessageBox.information(self,"No Results","Run docking first."); return
        path, _ = QFileDialog.getSaveFileName(self,"Save ZIP","docking_results.zip","ZIP (*.zip)")
        if not path: return
        cols = ['Ligand','Affinity (kcal/mol)','RMSD (Å)','Heavy Atoms','Ligand Efficiency']
        buf = io.BytesIO()
        with zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED) as zf:
            df = pd.DataFrame([{c: r.get(c,'') for c in cols} for r in good])
            zf.writestr('docking_results.csv', df.to_csv(index=False))
            if self._protein_data:
                zf.writestr(f'receptor.{self._protein_fmt}', self._protein_data)
            for r in good:
                safe = ''.join(c for c in r['Ligand'] if c.isalnum() or c in '._-').strip()
                zf.writestr(f'{safe}.pdbqt', r.get('PDBQT',''))
        with open(path,'wb') as fh: fh.write(buf.getvalue())
        QMessageBox.information(self,"Exported",f"✔ ZIP saved to:\n{path}")


# ══════════════════════════════════════════════════════════════════════════════
#  INTERACTION ANALYSIS WORKER
# ══════════════════════════════════════════════════════════════════════════════

# Interaction type colors & labels
INTERACTION_STYLES = {
    'HBDonor':      {'color': '#3d8ef0', 'label': 'H-Bond Donor',    'ls': '-',  'lw': 2.0},
    'HBAcceptor':   {'color': '#22c55e', 'label': 'H-Bond Acceptor', 'ls': '-',  'lw': 2.0},
    'Hydrophobic':  {'color': '#f59e0b', 'label': 'Hydrophobic',     'ls': '--', 'lw': 1.5},
    'PiStacking':   {'color': '#7c3aed', 'label': 'π–π Stacking',    'ls': '-',  'lw': 1.8},
    'PiCation':     {'color': '#ec4899', 'label': 'π–Cation',        'ls': '-',  'lw': 1.8},
    'CationPi':     {'color': '#ec4899', 'label': 'Cation–π',        'ls': '-',  'lw': 1.8},
    'Cationic':     {'color': '#ef4444', 'label': 'Ionic (+)',        'ls': '-',  'lw': 2.0},
    'Anionic':      {'color': '#ef4444', 'label': 'Ionic (−)',        'ls': '-',  'lw': 2.0},
    'VdWContact':   {'color': '#8494a8', 'label': 'VdW Contact',     'ls': ':',  'lw': 1.2},
    'XBAcceptor':   {'color': '#06b6d4', 'label': 'Halogen Bond',    'ls': '-',  'lw': 1.8},
}


class InteractionWorker(QThread):
    """Runs ProLIF protein-ligand interaction fingerprint in background."""
    finished = pyqtSignal(object)   # dict: {df, interactions, ligand_smiles}
    progress = pyqtSignal(str)
    error    = pyqtSignal(str)

    def __init__(self, protein_pdbqt: str, ligand_pdbqt: str, ligand_name: str):
        super().__init__()
        self.protein_pdbqt = protein_pdbqt
        self.ligand_pdbqt  = ligand_pdbqt
        self.ligand_name   = ligand_name

    def run(self):
        import warnings, tempfile, os, subprocess
        warnings.filterwarnings('ignore')
        try:
            import prolif
            import MDAnalysis as mda
            from rdkit import Chem
            from rdkit.Chem import AllChem

            self.progress.emit("[Interactions] Converting structures…")

            # Write protein PDBQT → temp file, convert to PDB
            with tempfile.NamedTemporaryFile(suffix='.pdbqt', mode='w', delete=False) as tf:
                tf.write(self.protein_pdbqt); prot_pdbqt_path = tf.name
            prot_pdb_path = prot_pdbqt_path.replace('.pdbqt', '_prot.pdb')
            subprocess.run([OBABEL, '-ipdbqt', prot_pdbqt_path, '-opdb',
                            '-O', prot_pdb_path], capture_output=True)

            # Write ligand PDBQT (first pose only) → convert to SDF
            with tempfile.NamedTemporaryFile(suffix='.pdbqt', mode='w', delete=False) as tf:
                # Take only first model (first pose)
                lines = []
                in_model = False
                for line in self.ligand_pdbqt.splitlines():
                    if line.startswith('MODEL'):
                        if in_model: break   # only first model
                        in_model = True
                    lines.append(line)
                    if line.startswith('ENDMDL'):
                        break
                tf.write('\n'.join(lines)); lig_pdbqt_path = tf.name

            lig_sdf_path = lig_pdbqt_path.replace('.pdbqt', '_lig.sdf')
            subprocess.run([OBABEL, '-ipdbqt', lig_pdbqt_path, '-osdf',
                            '-O', lig_sdf_path, '-h'], capture_output=True)

            if not os.path.exists(prot_pdb_path) or not os.path.exists(lig_sdf_path):
                self.error.emit("Structure conversion failed. Check OpenBabel installation.")
                return

            self.progress.emit("[Interactions] Running ProLIF fingerprint…")

            # Load with MDAnalysis + ProLIF
            u = mda.Universe(prot_pdb_path)
            lig_mol = Chem.SDMolSupplier(lig_sdf_path, removeHs=False)[0]
            if not lig_mol:
                self.error.emit("Could not parse ligand SDF."); return

            # Build prolif molecules
            prot_mol  = prolif.Molecule.from_mda(u)
            lig_plmol = prolif.Molecule.from_rdkit(lig_mol)

            # Run fingerprint
            fp = prolif.Fingerprint(count=False)
            fp.run_from_iterable([lig_plmol], prot_mol)

            df = fp.to_dataframe()

            # Parse interactions per residue
            interactions = {}  # {residue_str: [interaction_types]}
            if not df.empty:
                for col in df.columns:
                    val = df[col].iloc[0]
                    if val:
                        # col format: ('ligand', 'residue', 'interaction')
                        if len(col) == 3:
                            _, resname, itype = col
                            res_str = str(resname)
                            if res_str not in interactions:
                                interactions[res_str] = []
                            interactions[res_str].append(str(itype))

            # Get ligand SMILES for 2D drawing
            lig_smi = Chem.MolToSmiles(Chem.RemoveHs(lig_mol)) if lig_mol else ''

            # Cleanup
            for p in [prot_pdbqt_path, prot_pdb_path, lig_pdbqt_path, lig_sdf_path]:
                try: os.unlink(p)
                except: pass

            self.progress.emit(f"[Interactions] ✔ {len(interactions)} interacting residues found")
            self.finished.emit({
                'df':           df,
                'interactions': interactions,
                'ligand_smi':   lig_smi,
                'ligand_name':  self.ligand_name,
            })

        except Exception:
            self.error.emit(traceback.format_exc())


# ══════════════════════════════════════════════════════════════════════════════
#  2D INTERACTION DIAGRAM  (matplotlib radial network)
# ══════════════════════════════════════════════════════════════════════════════

class InteractionDiagram(FigureCanvas):
    """2D radial diagram showing protein residues and interaction types."""

    BG = C['bg']; AX = C['surface']; TXT = C['text']; TXT2 = C['text2']

    def __init__(self, parent=None):
        self.fig = Figure(facecolor=self.BG, figsize=(8, 8))
        super().__init__(self.fig)
        self.setMinimumSize(500, 500)
        self._placeholder()

    def _placeholder(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111, aspect='equal')
        ax.set_facecolor(self.BG)
        ax.text(0.5, 0.5, 'Run docking and select a result\nto view interactions',
                ha='center', va='center', color=self.TXT2,
                fontsize=12, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        self.fig.tight_layout(); self.draw()

    def render(self, interactions: dict, ligand_name: str, ligand_smi: str = ''):
        """
        interactions: {residue_str: [interaction_type, ...]}
        """
        import math
        from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
        from matplotlib.lines import Line2D
        import matplotlib.patheffects as pe

        self.fig.clear()
        ax = self.fig.add_subplot(111, aspect='equal')
        ax.set_facecolor(self.BG)
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.set_xticks([]); ax.set_yticks([])

        if not interactions:
            ax.text(0.5, 0.5, 'No interactions found in this pose.',
                    ha='center', va='center', color=self.TXT2,
                    fontsize=11, transform=ax.transAxes)
            self.fig.tight_layout(); self.draw(); return

        residues = list(interactions.keys())
        n = len(residues)
        R = 3.8   # radius of residue circle
        cx, cy = 0.0, 0.0

        # ── Ligand node (center) ─────────────────────────────────────────────
        # Try 2D structure from SMILES
        lig_drawn = False
        if ligand_smi:
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw, AllChem
                from io import BytesIO
                from matplotlib.image import imread
                import numpy as np

                mol2d = Chem.MolFromSmiles(ligand_smi)
                AllChem.Compute2DCoords(mol2d)
                img = Draw.MolToImage(mol2d, size=(160, 160),
                                      kekulize=True,
                                      options=Draw.MolDrawOptions())
                buf = BytesIO(); img.save(buf, format='PNG'); buf.seek(0)
                img_arr = imread(buf)
                # Place as image at center
                ax.imshow(img_arr, extent=[-1.2, 1.2, -1.2, 1.2],
                          zorder=5, aspect='auto')
                lig_drawn = True
            except Exception:
                pass

        if not lig_drawn:
            lig_circle = Circle((cx, cy), 0.9, facecolor=C['surface3'],
                                edgecolor=C['accent'], linewidth=2, zorder=5)
            ax.add_patch(lig_circle)
            ax.text(cx, cy, ligand_name[:12], ha='center', va='center',
                    color=C['accent'], fontsize=9, fontweight='bold', zorder=6)

        # ── Draw residue nodes and interaction edges ─────────────────────────
        for i, res in enumerate(residues):
            angle = 2 * math.pi * i / n - math.pi / 2
            rx = cx + R * math.cos(angle)
            ry = cy + R * math.sin(angle)

            itypes = interactions[res]

            # Determine dominant color (first interaction type)
            dom_itype = itypes[0] if itypes else 'VdWContact'
            style = INTERACTION_STYLES.get(dom_itype, INTERACTION_STYLES['VdWContact'])
            edge_color = style['color']

            # Draw edge(s) — one per interaction type
            for itype in itypes:
                st = INTERACTION_STYLES.get(itype, INTERACTION_STYLES['VdWContact'])
                # Direction from ligand center to residue, stopping before circles
                dx = rx - cx; dy = ry - cy
                dist = math.hypot(dx, dy)
                ux, uy = dx/dist, dy/dist
                start_r = 1.2; end_r = dist - 0.55
                x1 = cx + ux * start_r; y1 = cy + uy * start_r
                x2 = cx + ux * end_r;   y2 = cy + uy * end_r
                ax.plot([x1, x2], [y1, y2],
                        color=st['color'], lw=st['lw'],
                        ls=st['ls'], alpha=0.75, zorder=3)

            # Residue box
            box_w, box_h = 1.15, 0.45
            fancy = FancyBboxPatch(
                (rx - box_w/2, ry - box_h/2), box_w, box_h,
                boxstyle="round,pad=0.05",
                facecolor=C['surface3'], edgecolor=edge_color,
                linewidth=1.8, zorder=6
            )
            ax.add_patch(fancy)

            # Residue label
            # Parse residue string, e.g. "GLU128.A" → "GLU128"
            res_short = str(res).replace('.', '\n') if '.' in str(res) else str(res)
            ax.text(rx, ry, res_short, ha='center', va='center',
                    color=self.TXT, fontsize=7.5, fontweight='600',
                    zorder=7, linespacing=1.2)

            # Interaction type mini-dots along the edge (visual markers)
            mid_x = (x1 + x2) / 2; mid_y = (y1 + y2) / 2
            for k, itype in enumerate(itypes):
                st = INTERACTION_STYLES.get(itype, INTERACTION_STYLES['VdWContact'])
                offset = (k - len(itypes)/2 + 0.5) * 0.18
                px = mid_x - uy * offset; py = mid_y + ux * offset
                ax.plot(px, py, 'o', color=st['color'],
                        markersize=5, zorder=8, alpha=0.9)

        # ── Legend ───────────────────────────────────────────────────────────
        present_types = set()
        for ilist in interactions.values():
            present_types.update(ilist)

        legend_handles = []
        for itype in INTERACTION_STYLES:
            if itype in present_types:
                st = INTERACTION_STYLES[itype]
                legend_handles.append(
                    Line2D([0], [0], color=st['color'], lw=st['lw'],
                           ls=st['ls'], label=st['label'],
                           marker='o', markersize=5)
                )
        if legend_handles:
            leg = ax.legend(handles=legend_handles, loc='lower right',
                            fontsize=7, framealpha=0.85,
                            facecolor=C['surface2'], edgecolor=C['border2'],
                            labelcolor=self.TXT)

        ax.set_xlim(cx - R - 1.8, cx + R + 1.8)
        ax.set_ylim(cy - R - 1.8, cy + R + 1.8)
        ax.set_title(f"Protein–Ligand Interactions: {ligand_name}",
                     color=self.TXT, fontsize=11, fontweight='600', pad=10)
        self.fig.tight_layout()
        self.draw()


# ── Entry Point ────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import sys

    def _eh(t, v, tb):
        if issubclass(t, KeyboardInterrupt): sys.__excepthook__(t,v,tb); return
        try: QMessageBox.critical(None, f"Error — {t.__name__}", ''.join(__import__('traceback').format_exception(t,v,tb))[:2000])
        except: sys.__excepthook__(t,v,tb)
    sys.excepthook = _eh

    app = QApplication(sys.argv)
    app.setStyleSheet(QSS)
    from PyQt6.QtGui import QPalette, QColor as QC2
    pal = QPalette()
    pal.setColor(QPalette.ColorRole.Window,          QC2(C['bg']))
    pal.setColor(QPalette.ColorRole.WindowText,      QC2(C['text']))
    pal.setColor(QPalette.ColorRole.Base,            QC2(C['surface2']))
    pal.setColor(QPalette.ColorRole.AlternateBase,   QC2(C['surface']))
    pal.setColor(QPalette.ColorRole.Text,            QC2(C['text']))
    pal.setColor(QPalette.ColorRole.Button,          QC2(C['surface3']))
    pal.setColor(QPalette.ColorRole.ButtonText,      QC2(C['text']))
    pal.setColor(QPalette.ColorRole.Highlight,       QC2(C['accent_d']))
    pal.setColor(QPalette.ColorRole.HighlightedText, QC2('#f0f4f9'))
    app.setPalette(pal)

    win = EasyDockingApp()
    win.show()
    sys.exit(app.exec())
