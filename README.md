# rpGlobalScore

Calculate the global score based on the thermodynamics, FBA, Selenzyme analysis. For more information of the way the global score is calculated please refer to the following document TODO. The tool returns the top scoring pathways as set by the user. 

## Information Flow

### Input

Required information:
* Single rpSBML or collection of SBML file as a tar.xz
* Maximal number of steps (as run in RetroPath2.0)
* Number of steps weight
* Selenzyme weight
* FBA weight
* Thermodynamics weight
* Number of best scoring pathways to return

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

* Docker - [Install](https://docs.docker.com/v17.09/engine/installation/)
* libSBML - [Anaconda library](https://anaconda.org/SBMLTeam/python-libsbml)

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Version 0.1

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson
