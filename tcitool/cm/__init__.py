import glob, re, os
import xml.etree.ElementTree
import numpy as np
import matplotlib

def _process_color(color_str,alpha=None):
    _rbg_match = re.match(r'^rgb\(([0-9.]+),([0-9.]+),([0-9.]+)\)$',color_str)
    _rgb_val = np.array([np.nan,np.nan,np.nan])
    if alpha is not None:
        _rgb_val = np.append(_rgb_val,np.clip(float(alpha),0.,1.))
    if _rbg_match is not None:
        for _i in range(3):
            _group = _rbg_match.groups()[_i]
            if '%' in _group:
                _v = np.clip(int(_group[:-1]),0,100)/100
            else:
                _v = np.clip(int(_group),0,255)/255
            _rgb_val[_i] = _v
        return _rgb_val
    _hex3_match = re.match(r'^\#[0-9A-Fa-f]{3}$',color_str)
    if _hex3_match is not None:
        _,_r,_g,_b = color_str
        color_str = '#'+_r+_r+_g+_g+_b+_b
    _hex6_match = re.match(r'^\#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$',color_str)
    if _hex6_match is not None:
        _rgbhex = _hex6_match.groups()
        _rgbdec = np.array(list(map(lambda s: int(s.lower(), 16), _rgbhex)))
        _rgb_val[0],_rgb_val[1],_rgb_val[2] = _rgbdec[0]/255,_rgbdec[1]/255,_rgbdec[2]/255
        return _rgb_val
    raise ValueError('Color format not recognised')

try:
    cmapnames = []
    cmap_d = {}
    _xml_namespaces = {'svg':'http://www.w3.org/2000/svg'}
    folder = os.path.dirname(os.path.realpath(__file__))
    for _file in glob.iglob(os.path.join(folder,'*.svg')):
        _tree = xml.etree.ElementTree.parse(_file)
        for _gradient in _tree.findall('.//svg:linearGradient',_xml_namespaces):
            _cmap_name = _gradient.attrib["id"]
            _cmap_name_r = _cmap_name+'_r'
            _cmap_array = []
            _cmap_offset = []
            for _step in _gradient:
                _step_offset, _step_color, _step_alpha = (_step.attrib[p] for p in ['offset','stop-color','stop-opacity'])
                _cmap_array.append(_process_color(_step_color,_step_alpha))
                _cmap_offset.append(np.round(float(_step_offset[:-1])/100,5))
            _cmap_offset[0],_cmap_offset[-1] = 0.,1.
            _cmap = matplotlib.colors.LinearSegmentedColormap.from_list(name=_cmap_name,colors=list(zip(_cmap_offset,_cmap_array)))
            _cmap_r = _cmap.reversed(name=_cmap_name_r)
            cmapnames.append(_cmap_name)
            cmap_d[_cmap_name] = _cmap
            cmap_d[_cmap_name_r] = _cmap_r
            locals()[_cmap_name] = _cmap
            locals()[_cmap_name_r] = _cmap_r
    del _xml_namespaces, _file, _tree, _gradient, _cmap_name, _cmap_array, _cmap_offset, _step_offset, _step_color, _step_alpha, glob, re, os, xml, np, matplotlib
except:
    raise