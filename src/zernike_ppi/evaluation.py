import numpy as np
from scipy import metrics
def evaluate(scores, labels, k=20):
 scores=np.asarray(scores,float); labels=np.asarray(labels,int); order=np.argsort(scores); y=labels[order];
 return {"ROC_AUC":float(metrics.roc_auc_score(labels,-scores)) if len(np.unique(labels))>1 else float("nan"),"PR_AUC":float(metrics.average_precision_score(labels,-scores)) if labels.any() else float("nan"),"precision_at_k":float(y[:k].mean()) if len(y) else float("nan"),"recall_at_k":float(y[:k].sum()/max(1,y.sum())) if len(y) else float("nan"),"mean_rank":float(np.mean(np.where(labels[order]>0)[0]+1)) if labels.any() else float("nan")}
