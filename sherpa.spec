# -*- mode: python -*-

block_cipher = None


a = Analysis(['sherpa.py'],
             pathex=['D:\\WORK\\projects\\1_urbIam\\1_CODE_MATLAB\\SHERPA\\SHERPA-GITHUB\\Sherpa'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='sherpa',
          debug=False,
          strip=False,
          upx=False,
          runtime_tmpdir=None,
          console=True )
