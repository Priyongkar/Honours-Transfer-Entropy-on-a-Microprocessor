# Sim2600, ported to Python 3, with clamps

Upstream: https://github.com/ericmjonas/Sim2600 (CC0), commit f7198770563073f8c38cb54f3937f091cbeb503f.
The ROMs are in the upstream repository and are not copied here.

The patch:
1. ports the code to Python 3 (lib2to3; pickle loaded with `encoding="latin1"`; relative imports;
   ROM bytes read as integers);
2. builds the Cython engine with `setup_py3.py` (C++11, language level 3);
3. adds forced wires: `_forced[w] = 1` (low) or `2` (high) makes the wire's group take that value during
   settling. Unforced behaviour is unchanged: the regenerated Donkey Kong recording equals the published
   file in every analysed sample.

Apply: `git clone https://github.com/ericmjonas/Sim2600 && cd Sim2600 && git checkout f719877 && git apply ../sim2600_python3_port_and_forced_wire_clamps.patch && cp ../setup_py3.py . && python setup_py3.py build_ext --inplace`
