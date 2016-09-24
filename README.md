this is the implementation of ptychography algorithm, a very simple toy program

it will make two dir in your present work dir: det_mod and recon
det_mod is the modulus of diffraction pattern get by the detector at
different position.
recon is the reconstruct image generate during the iteration.

startPos and endPos, this two args will control the overlap rate

m is the radius of support

You can replace the HIO2D_ver1 with ER2D_ver1 or other phase retrieval algorithm.

mail: csu_scj2015@foxmail.com