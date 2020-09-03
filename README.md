# Medchem design moves for generating new compounds

This repository contains a python library for exploring the medicinal chemistry design moves, as reported in publication: "The playbook of Medicinal Chemistry Design Moves, xxx-xxx-xxx"
The library is adpated from original [mmpdb code](https://github.com/rdkit/mmpdb)

## Dependency

  - python >=3
  - rdkit

## Dataset

  - please download the design moves (transformation) dataset from [here](www.google.com). We provide the design moves from ChEMBL database at radius 3.

## How to run programm
```sh
python mmpdb mmpCompoundGenerator  --tsmiles "O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N" --transformation_db "chemblDB3.sqlitdb" --replaceGroup "*S(=O)(=O)(N)" --tradius 3 --toutput output.txt --tmin-pairs 100
```

## Output

|original_smi|transformed_smi|original_frag|new_frag|rule_freq|ex_lhs_cpd_id|ex_rhs_cpd_id|
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccccc2)cc1|[*]S(N)(=O)=O|[*:1][H]|1103|CHEMBL51385|CHEMBL1201104|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|COc1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]OC|607|CHEMBL51385|CHEMBL8441|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(Cl)cc2)cc1|[*]S(N)(=O)=O|[*:1]Cl|523|CHEMBL51385|CHEMBL462|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(C)(=O)=O)cc2)cc1|[*]S(N)(=O)=O|[*:1]S(C)(=O)=O|497|CHEMBL468367|CHEMBL507789|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(F)cc2)cc1|[*]S(N)(=O)=O|[*:1]F|491|CHEMBL3426428|CHEMBL3426432|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]C|408|CHEMBL51385|CHEMBL274877|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(C(=O)O)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(=O)O|275|CHEMBL1966874|CHEMBL2094690|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc([N+](=O)[O-])cc2)cc1|[*]S(N)(=O)=O|[*:1][N+](=O)[O-]|265|CHEMBL51385|CHEMBL8682|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(O)cc2)cc1|[*]S(N)(=O)=O|[*:1]O|255|CHEMBL3632832|CHEMBL1341020|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(C(N)=O)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(N)=O|215|CHEMBL3901141|CHEMBL3973258|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(Br)cc2)cc1|[*]S(N)(=O)=O|[*:1]Br|205|CHEMBL51385|CHEMBL450762|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(C#N)cc2)cc1|[*]S(N)(=O)=O|[*:1]C#N|193|CHEMBL3426428|CHEMBL3426430|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|CC(=O)Nc1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]NC(C)=O|170|CHEMBL1490019|CHEMBL1716793|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|CC(=O)c1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(C)=O|152|CHEMBL2163818|CHEMBL1206418|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|CCOc1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]OCC|141|CHEMBL2163818|CHEMBL1206420|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(C(F)(F)F)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(F)(F)F|137|CHEMBL51385|CHEMBL501107|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|COC(=O)c1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(=O)OC|136|CHEMBL3695775|CHEMBL3695774|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|CCOC(=O)c1ccc(-n2nc(C(F)(F)F)cc2-c2ccc(C)cc2)cc1|[*]S(N)(=O)=O|[*:1]C(=O)OCC|132|CHEMBL285831|CHEMBL241971|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(N)cc2)cc1|[*]S(N)(=O)=O|[*:1]N|126|CHEMBL51385|CHEMBL463|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(N3CCOCC3)cc2)cc1|[*]S(N)(=O)=O|[*:1]N1CCOCC1|114|CHEMBL3679534|CHEMBL3679542|
|O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N|Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(N(C)C)cc2)cc1|[*]S(N)(=O)=O|[*:1]N(C)C|108|CHEMBL3695775|CHEMBL3695770|


* note that output file contains more columns than shown above. For clearity only main columns are shown.

## Input Options 

|parameter|meaning |
|---------|--------|
|--transformation_db| transformation database file|
|--tsmiles| input query molecule (in smiles format)|
|--replaceGroup|fragment to be replace in a query molecule (in smiles format). Note that it must contains * or [*] to indicate attachment point. 1) single cut example *N1CCOCC1, 2) double cut example *N1CCOC(*)C1
|--tradius|chemical envionment radius for design move. It must equal to radius of the transformation database. So in this case 3.|
|--toutput|name of output text file|
|--tmin-pairs|consider design moves that have at-least tmin-pairs examples. This is in other words the freqency of a design move. For instacne if we set tmin-pairs to 5: it say that consider all design moves that are derived from at-least five MMP pairs.| 


## Outpt Explnation
|column|meaning |
|---------|--------|
|original_smi|smiles of a query molecule|
|transformed_smi|a new design molecule|
|constant_smi|smiles of constant part i.e. part which did not change|
|original_frag|fragment in a query to replace|
|new_frag|new suggested fragment|
|envsmi|chemical environment of replacementi (design move)|
|rule_freq|frequency of design move (transformation)|
|ex_lhs_cpd_id|an example of MMP (left hand compound id) from ChEMBL|
|ex_rhs_cpd_id|an example of MMP (right hand compound id) from ChEMBL|
