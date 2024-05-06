Neutrino flavor evolution with advection: simulating dense astrophysical environments in 2D
====

Description: We solve the schõdinger equation in 2D to simulate how a dense ensemble of neutrinos oscillates in dense astrophysical environments. 
These results are relevant for core-collapse supernovae and the mergers of neutron stars

## Journal reference: 
[Shashank Shalgar et al JCAP06(2020)048](https://iopscience.iop.org/article/10.1088/1475-7516/2020/06/048), [ePrint:1911.09110](https://arxiv.org/abs/1911.09110)

## Figures & animations that appear in my publication:

![imagesums](https://github.com/ianpaga/neutrino_gas_advection/assets/57350668/64e0b0c3-ae45-4a34-a17a-f7261e1facef)
![image001](https://github.com/ianpaga/neutrino_gas_advection/assets/57350668/564dc7d7-9c34-4fe7-ac34-caab4468751a)
![image100](https://github.com/ianpaga/neutrino_gas_advection/assets/57350668/d7a83650-75b1-47c7-b828-ad07b93f09de)
![image200](https://github.com/ianpaga/neutrino_gas_advection/assets/57350668/9628ecbc-e643-4d73-8bb6-efaef060e4ae)

## Requirements:
- C++ compiler
- [Boost Library](https://www.boost.org/)
- OpenMP
- Python, Matplotlib, NumPy

## Compiling and running:
- Run ./chpc to compile
- Run executable *.out
- Outputs *.raw and *.sum are large files. Use the bash script splitandplot.sh to slice the data into smaller files and make plots.
- The splitandplot.sh from the previous step  calls the Python script plotaframe.py
