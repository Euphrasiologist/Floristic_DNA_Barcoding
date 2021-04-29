# Dating the British flora phylogeny

I'm looking to do this to make the DNA Barcoding paper Alex and I have been working on more robust. I think the idea came from Pete Hollingsworth. It's a good idea to do this anyway, because we will need an accurate, dated phylogeny for a potential project on chromosome number and genome size evolution incorporating this data.

## Making the tree

Alex suggested https://github.com/blackrim/treePL to date the tree(s). He also suggested https://doi.org/10.1038/s41559-020-1241-3 to poach the fossil calibrations. 

To run:

`treePL config.treepl`

Some useful documentation: https://arxiv.org/pdf/2008.07054.pdf

## Calibration points

Let's just put this straight into CSV format in case. There's a max of two rows per clade, as these are the min/max values.

Clade, MYA, Fossil
Angiospermae, 129.4, †Retisulc-Muriverm
Angiospermae, 132.1, †earliest angiosperm-like pollen grains
Monocotyledoneae, 113, †Liliacidites sp. A
Eudicotelydonae, 125, †Hyrcantha decussata
Nymphaeaceae, 110.8, †Monetianthus mirus
Ceratophyllaceae, 100.5, †Donlesia dakotensis
Ceratophyllaceae, 127.2, †Montsechia vidalii
Alismataceae, 23.03, †Caldesia brandoniana
Alismataceae, 72.1, †Cardstonia toldmanii
Araceae, 47, †Araciphyllites tertiarius
Araceae, 123, †Araciphyllites tertiarius
Orchidaceae, 23.2, †Dendrobium winikaphyllum
Orchidaceae, 15, †Meliorchis caribea
Poaceae, 30.44, †Leersia seifhennersdorfensis
Poaceae, 66, †Changii indicum
Ranunculaceae, 56, †Paleoactaea nagelii
Malvaceae, 33.9, †Florissantia ashwillii
Malvaceae, 72.1, †Javelinoxylon weberi
Fabaceae, 15, †Acacia eocaribbeanensis
Fabaceae, 56, †Diplotropis-like leaves and fruits
Salicaceae, 48.5, †Populus tidwelli
Polygonaceae, 66, †Polygonocarpum johnsonii
Ericaceae, 66, †Leucothoe praecox
Apiaceae, 66, †Carpites ulmiformis
Apiaceae, 41.2, †Umbelliferospermum latahense
Asteraceae, 47.46, †Raiguenrayun cura
Asteraceae, 72.1, †Tubulifloridites lilliei type A
Lamiaceae, 27.82, †Ajuginucula smithii
Solanaceae, 52.22, †Physalis infinemundi
Solanaceae, 33.9, †Solanispermum reniforme