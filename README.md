
In the folder “feature test” are the features of the residue pairs on tetramers in testing set. Each tetraner has four chains, which means that there are six possible interfaces. For example, protein 1bv4 has four chains: A, B, C and D. In file 1bv4_AB.txt are features of residue pairs from chain A and B of protein 1bv4.


The features of only one tetramer are given here. Full test set data can be obtained from this link: [ftp://202.112.126.135/pub/Tetramer/](ftp://202.112.126.135/pub/Tetramer/)


## Requirements

- Python 3

- numpy\==1.16.0

- pandas==1.0.1

- tensorflow==1.14.0

- scipy==1.4.1

- h5py==2.10.0

  You can install packages by :
  ```
  pip install -r requirements.txt
  ```
  
  

## Run

Take 1bv4_AB.txt as an example:
```
python model_test.py 1bv4_AB.txt
```
Then you can get result of 1bv4_AB in folder “results”.

The first two columns of the result are the corresponding sequence numbers for each residue pair. The third column is the label of each residue pair. “1.0” means it is an interface residue pair, “0.0” means not. 

The order of the results from top to bottom is arranged according to the predicted probability on the interface from large to small. The more occurrences of “1” in the top result sequence represents the better results. 









  
