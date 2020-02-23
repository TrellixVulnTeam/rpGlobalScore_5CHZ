# rpGlobalScore

Calculate the global score based on the thermodynamics, FBA, Selenzyme analysis. For more information of the way the global score is calculated please refer to the following document TODO. The tool returns the top scoring pathways as set by the user. 

## Information Flow

### Input

Required information:
* **Input rpSBML**: Single rpSBML or collection of SBML file as a tar.xz
* **Number of best scoring pathways to return**

Advanced Options:
* **Name of the heterologous pathway**: (default: rp_pathway) Groups ID of the heterologous pathway
* **Thermodynamics ID**: (default: dfG_prime_m)
* **Thermodynamics normalisation ceiling**: (default: 8901.2)
* **Thermodynamics normalisation floor**: (default: -7570.2)
* **FBA objective ID**: (default: obj_RP1_sink__restricted_biomass)
* **FBA normalisation ceiling**: (default: 3.0)
* **FBA normalisation floor**: (default: 0.0)
* **IP address of the rpGlobalScore REST service**: IP address of the REST service
* **Maximal number of steps (as run in RetroPath2.0)**
* **Number of steps weight**
* **Selenzyme weight**
* **FBA weight**
* **Thermodynamics weight**

### Output

* **rpGlobalScore**: Single or collection of rpSBML

## Installing

To compile the docker use the following command:

```
docker build -t brsynth/rpglobalscore-rest:dev -f Dockerfile .
```

And then run the container with the follwing command:

```
docker run -p 8881:8888 brsynth/rpglobalscore-rest:dev
```

### Prerequisites

=======

Required information:
* **Input/Output**: Single rpSBML or collection of SBML file as a tar.xz
* **Maximal number of steps (as run in RetroPath2.0)**
* **Number of steps weight**
* **Selenzyme weight**
* **FBA weight**
* **Thermodynamics weight**
* **Number of best scoring pathways to return**

Advanced Options:
* **Name of the heterologous pathway**: (default: rp_pathway) Groups ID of the heterologous pathway
* **Thermodynamics ID**: (default: dfG_prime_m)
* **Thermodynamics normalisation ceiling**: (default: 8901.2)
* **Thermodynamics normalisation floor**: (default: -7570.2)
* **FBA objective ID**: (default: obj_RP1_sink__restricted_biomass)
* **FBA normalisation ceiling**: (default: 3.0)
* **FBA normalisation floor**: (default: 0.0)
* **IP address of the rpGlobalScore REST service**: IP address of the REST service

### Output

* **rpGlobalScore**: Collection (as tar.xz) or single rpSBML file

## Installing

To compile the docker use the following command:

```
docker build -t brsynth/rpglobalscore-rest:dev -f Dockerfile .
```

And then run the container with the follwing command:

```
docker run -p 8881:8888 brsynth/rpglobalscore-rest:dev
```

## Algorithm

The following features are normalised: Length of the pathway, Gibbs free energy, FBA of the target (with fixed fraction of biomass) and Selenzyme scores.

### Min-Max Feature Scaling

![img](http://latex.codecogs.com/svg.latex?x%27%3D%5Cfrac%7Bx-%5Clfloor%08x%5Crfloor%7D%7B%5Clceil%08x%5Crceil-%5Clfloor%08x%5Crfloor%7D)

### Global Score

![img](http://latex.codecogs.com/svg.latex?%5Cbar%7Bx%7D%3D%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5En%7Bw_ix_i%7D%7D%7B%5Csum_%7Bi%3D1%7D%5En%7Bw_i%7D%7D)

## Prerequisites

* Docker - [Install](https://docs.docker.com/v17.09/engine/installation/)
* libSBML - [Anaconda library](https://anaconda.org/SBMLTeam/python-libsbml)

## Contributing

TODO

## Versioning

Version 0.1

## Authors

* **Melchior du Lac** 

### Acknowledgments

* Thomas Duigou
* Joan Hérisson

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
