# Predicting Similarity to Closest Known Reactions as a Surrogate for Reactivity

This repository contains the code and resources for **Predicting Similarity to Closest Known Reactions as a Surrogate for Reactivity**. The project aims to efficiently evaluate reaction feasibility by leveraging similarity-based metrics and machine learning models.

---

## 📦 1. Requirements
Make sure you have the following packages installed:
## 1.1 Package
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

## 🚀 2. How to Use
### 2.1 Reaction generation
You can find reaction generation part in folder **(generate_model)**. We provided a well trained model **uspto_model_v1_epoch60.pth**, so you can use it dirctly. We also provided the training dataset **uspto_reaction_smile_aug.txt**, you can also train a new model.
* genarationv6_mult_step.ipynb: it generated new molecules.
* combine_result_v1.ipynb: it combined new molecules to new reaction based on the reaction template.

### 2.2 Correlation evaluation
The Correlation Evaluation folder contains a comprehensive evaluation of the relationship between various reaction similarity measures and reactivity scores. Specifically, we assess the correlation of four different reaction similarity metrics and two reactivity scores provided by ASKCOS mentioned in the paper with the reactivity of 180 manually labeled Suzuki reactions.

### 2.3 Similarity prediction
You can find the XGBoost model to predict similarity and the SHAP analysis part. We apply SHAP (SHapley Additive exPlanations) to interpret the contribution of each feature in the similarity prediction process. This helps in understanding the model's decision-making and provides transparency in reaction feasibility prediction.

## 📄 3. License
This project is licensed under the MIT License.


