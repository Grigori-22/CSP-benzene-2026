Aby odtworzyć eksperyment należy:

1) Mając w folderze plik benzene.xyz użyć kodu PairGen.py. Wygeneruje on ponad 700 plików z rozszerzeniem .xyz zawierającymi struktury dimerów.

2) Pliki ze strukturami dimerów wykorzystać do obliczeń kwantowo-chemicznych w programie DFTB+. Przykładowe slurmy do takich obliczeń znajdują się w folderze pairs_slurms. Za pomocą np. "grep" uzyskać energie oddziaływań i zapisać je posortowane w pliku results.dat.

3) Mając pliki xyz z pierwszego skryptu oraz plik results.dat użyć kodu CrystalGen. On wygeneruje 1133 plików z rozszerzeniem .gen zawierające struktury krystaliczne.

4) Pliki z folderu crystals wykorzystać do obliczeń analizy wibracyjnej. Najpierw uzyskać Hessian, a następnie obliczyć drgania normalne. Przykładowe pliki do takich obliczeń w DFTB+ znajdują się folderze crystal_slurms
