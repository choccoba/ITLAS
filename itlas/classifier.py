"""Phase Classifier Module"""

import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier


class PhaseClassifier:
    """HBV infection phase classifier."""
    
    def __init__(self, model_type='random_forest'):
        self.model_type = model_type
        self.model = None
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
    
    def _get_model(self):
        if self.model_type == 'xgboost':
            try:
                from xgboost import XGBClassifier
                return XGBClassifier(n_estimators=100, max_depth=5, random_state=42)
            except ImportError:
                self.model_type = 'random_forest'
        return RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
    
    def _extract_features(self, adata):
        from .signature_score import calculate_it_signature_score
        from .fm_gsea import calculate_pathway_scores
        
        sig_scores = calculate_it_signature_score(adata, return_components=True)
        pathway_df = calculate_pathway_scores(adata)
        
        features = [np.array(list(sig_scores.values())).T, pathway_df.values]
        return np.hstack(features)
    
    def fit(self, adata, stage_col='Stage'):
        X = self._extract_features(adata)
        y = self.label_encoder.fit_transform(adata.obs[stage_col])
        X_scaled = self.scaler.fit_transform(X)
        self.model = self._get_model()
        self.model.fit(X_scaled, y)
        return self
    
    def predict(self, adata):
        X = self._extract_features(adata)
        X_scaled = self.scaler.transform(X)
        y_pred = self.model.predict(X_scaled)
        return self.label_encoder.inverse_transform(y_pred)
    
    def cross_validate(self, adata, stage_col='Stage', n_folds=5):
        X = self._extract_features(adata)
        y = self.label_encoder.fit_transform(adata.obs[stage_col])
        X_scaled = self.scaler.fit_transform(X)
        
        model = self._get_model()
        cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
        acc = cross_val_score(model, X_scaled, y, cv=cv, scoring='accuracy')
        f1 = cross_val_score(model, X_scaled, y, cv=cv, scoring='f1_weighted')
        
        return {'accuracy': f"{acc.mean():.3f}±{acc.std():.3f}",
                'f1_score': f"{f1.mean():.3f}±{f1.std():.3f}"}
