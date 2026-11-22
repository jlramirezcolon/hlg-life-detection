"""Make matplotlib use Helvetica Neue Light for all text and Helvetica Neue Bold for weight='bold' (panel letters)."""
from pathlib import Path
import matplotlib
from matplotlib import font_manager as fm

HERE = Path(__file__).resolve().parent


def use_helvetica_neue_light():
    fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if f.name != 'Helvetica Neue']
    for fn, w in (('HelveticaNeue-Light.ttf', 300), ('HelveticaNeue-Bold.ttf', 700)):
        if not (HERE / fn).exists():   # the fonts are not redistributed; on macOS, extract them from the system font
            import subprocess, sys
            face = {'HelveticaNeue-Light.ttf': 'Helvetica Neue Light', 'HelveticaNeue-Bold.ttf': 'Helvetica Neue Bold'}[fn]
            subprocess.run([sys.executable, str(HERE / 'extract_ttc_face.py'), '/System/Library/Fonts/HelveticaNeue.ttc',
                            face, str(HERE / fn)], check=True)
        fm.fontManager.addfont(str(HERE / fn))
        fm.fontManager.ttflist[-1].weight = w
    matplotlib.rcParams.update({'font.family': 'Helvetica Neue', 'font.weight': 300, 'axes.titleweight': 300,
                                'axes.labelweight': 300, 'pdf.fonttype': 42, 'mathtext.fontset': 'custom',
                                'mathtext.rm': 'Helvetica Neue:light', 'mathtext.it': 'Helvetica Neue:light'})
    fm._get_font.cache_clear() if hasattr(fm, '_get_font') else None
