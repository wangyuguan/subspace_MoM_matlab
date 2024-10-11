# Subspace Method of Moments for Cryo-EM

A Complementary software for the paper: "Subspace Method of Moments for  Ab Initio 3-D Single-particle Cryo-EM Reconstruction”.

Paper authors: Jeremy Hoskins, Yuehaw Khoo, Oscar Mickelin, Amit Singer and Yuguan Wang.

Version: 1.0.

Release date: 8 Oct 2024.

Matlab code.

### Usage:

A simple example that can be completed in a few hours is available via "simple_example/run_test.m", after installing the following third party packages and downloading the following third party codes and data to "thirdparty/" folder.  

For a reproduction of the figures in the paper, see "plot_against_L/" and "plot_N_SNRs/" folders.

### 3rd party package, code and data:

This package uses the following softwares produced by 3rd parties:
1.	ASPIRE package, see https://github.com/PrincetonUniversity/aspire. 
2.	Manopt, see https://www.manopt.org.  

This package uses the following codes and data produced by 3rd parties:
1.	besselzero.m,		  see	[link](https://github.com/jinwar/matnoise/blob/master/besselzero.m).

2.	harmonicY.m,		  see [link](https://github.com/jmontalt/harmonicY/blob/master/harmonicY.m).

3.	lgwt.m,		  see [link](https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes).

4.	ReadMRC.m,		  see [link](https://github.com/nogaleslab/FreeHand/blob/master/ReadMRC.m).

5.	spherequad.m,		  see [link](https://people.sc.fsu.edu/~jburkardt/m_src/sphere_quad/sphere_quad.html).

6.	wignerD.m,		  see [link](https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fe5a37c32-4a80-11e4-9553-005056977bd0%2Fdea46a4f-38b6-68b6-2990-f52999540413%2Ffiles%2FwignerD.m&embed=web).

7.	WriteMRC.m,		  see [link](https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fch.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fe56f34df-4a80-11e4-9553-005056977bd0%2Fe5f1844b-3976-4384-be1d-b237e1f44f1c%2Ffiles%2FEMIODist2%2FWriteMRC.m&embed=web).

8.	WriteMRCHeader.m,		  see [link](https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fnl.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F50091%2Fversions%2F3%2Fcontents%2FUtils%2FwriteMRCHeader.m&embed=web).

9.	precomputed quadrature rules,		  see [link](https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/index.php.en).

10.	EMDB Data,    see [link](https://www.ebi.ac.uk/emdb/EMD-34948). 
