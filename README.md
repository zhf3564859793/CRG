# CRG

This is the code for 'Chemical Reaction Reactivity Prediction Based on Similarity '  
## 1. Package
* python==3.9.18  
* numpy==1.26.2  
* pandas==2.1.1 
* scikit-learn==1.3.2  
* rdkit==2023.9.2
* scipy==1.11.4
* seaborn==0.13.1
* tqdm==4.66.1
* xgboost==2.0.3

We suggest that you can install these packages by pip.

## 2. How to use
### 2.1 Reaction generation
You can find reaction generation part in folder(generate_model). We provided a well trained model uspto_model_v1_epoch60.pth, so you can use it dirctly(we also provided the training dataset “uspto_reaction_smile_aug.txt“,yo can also train a new model).
* genarationv6_mult_step.ipynb: it generated new molecules.
* combine_result_v1.ipynb: it combined new molecules to new reaction based on the reaction template.

### 2.2 Correlation evaluation
You can find different methods for correlation evaluation.

### 2.3 Similarity prediction
You can find the XGBoost model to predict similarity and the SHAP analysis part.
