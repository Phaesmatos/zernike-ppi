import numpy as np
def _auc(y, score):
    order=np.argsort(-score); y=np.asarray(y)[order]; pos=y.sum(); neg=len(y)-pos
    if not pos or not neg: return float('nan')
    ranks=np.where(y==1)[0]+1
    return float(1-(ranks.sum()-pos*(pos+1)/2)/(pos*neg))
def evaluate(scores, labels, k=20):
 scores=np.asarray(scores,float); labels=np.asarray(labels,int); order=np.argsort(scores); y=labels[order];
 pr=np.cumsum(y)/np.arange(1,len(y)+1); return {"ROC_AUC":_auc(labels,-scores) if len(np.unique(labels))>1 else float("nan"),"PR_AUC":float(np.sum(pr[y==1])/max(1,y.sum())) if labels.any() else float("nan"),"precision_at_k":float(y[:k].mean()) if len(y) else float("nan"),"recall_at_k":float(y[:k].sum()/max(1,y.sum())) if len(y) else float("nan"),"mean_rank":float(np.mean(np.where(labels[order]>0)[0]+1)) if labels.any() else float("nan")}
