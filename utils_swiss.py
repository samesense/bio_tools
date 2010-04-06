import PyMozilla

def getSwissModelPage(sequence):
    percent = 'missing'
    pdb_page = ''

    moz_emu = PyMozilla.MozillaEmulator(cacher=None, trycount=0)
    page = moz_emu.download('http://swissmodel.expasy.org/repository/')

    fields = [['query_1_input', sequence],
              ['query_1_type', 'SwissProt/TrEMBL AC'],
              ['db_sprot', 'swissprot'],
              ['query_limit', '50'],
              ['Search', 'SEARCH']]
    try:
        page = moz_emu.post_multipart('http://swissmodel.expasy.org/repository/?pid=smr03&uid=&token=&zid=async',
                                      fields, [])        
        this_line_is_percent = False
        suffix = ''
        
        for line in page.split('\n'):
            if this_line_is_percent:
                percent = line.split('>')[1].split('%')[0]
                this_line_is_percent = False
            elif line.find('Sequence identity') != -1:
                this_line_is_percent = True
            elif line.find('download') != -1:
                if line.find('PDB') != -1:
                    suffix = line.split('?pid')[1].split("\"")[0]
                    break
    
        if suffix != '':
            pdb_page = moz_emu.download('http://swissmodel.expasy.org/repository/?pid'
                                        + suffix)
    except: pass

    return [percent, pdb_page]

