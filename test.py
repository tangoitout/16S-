from qiime2 import Artifact
import pandas as pd

df=Artifact.load("percentile_normalized.qza").view(pd.DataFrame)
converted=Artifact.import_data("FeatureTable[Frequency]",df)
converted.save("pnorm_freq.qza")