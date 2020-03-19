# rpGlobalScore

Tool that reads a collection of rpSBML files (in a tar.xz) and calculates a global score based on a series of weights given by the user

## Getting Started

This is a docker galaxy tools, and thus, the docker needs to be built locally where Galaxy is installed. 

## Input

Required information:
* **-input**: (string) Path to either tar.xz input collection of rpSBML files or a single rpSBML file.
* **-input_format**: (string) Format of the input

Advanced options:
* **-fba_ceil**: (float, default=3.0) FBA ceiling
* **-fba_floor**: (float, default=0.0) FBA floor
* **-thermo_ceil**: (float, default=8901.2) Thermodynamics ceiling
* **-thermo_floor**: (float, default=-7570.2) Thermodynamics floor
* **-weight_rp_steps**: (float, default=0.0) Number of steps weight
* **-max_rp_steps**: (integer, default=15) Maximal number of steps (as run in RetroPath2.0)
* **-weight_rule_score**: (float, default=0.0) Reation rule score weight
* **-weight_fba**: (float, default=0.699707) FBA weight
* **-weight_thermo**: (float, default=0.8334961) Pathway thermodynamics weight
* **-pathway_id**: (string, default=rp_pathway) Name of the heterologous pathway
* **-objective_id**: (string, default=obj_RP1_sink__restricted_biomass) Name of the heterologous pathway objective function
* **-thermo_id**: (string, default=dfG_prime_m) Name of the thermodynamics

## Output

* **output**: (string) Path to the output file

## Dependencies

* Base Docker Image: [brsynth/rpbase](https://hub.docker.com/r/brsynth/rpbase)

## Installing

To build the image using the Dockerfile, use the following command:

```
docker build -t brsynth/rpglobalscore-standalone:dev -f Dockerfile .
```

### Running the tests

To run the test, untar the test.tar.xz file and run the following command:

```
python run.py -input test/test_rpGlobalScore.tar -input_format tar -output test/test_output.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

v0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpGlobalScore?
