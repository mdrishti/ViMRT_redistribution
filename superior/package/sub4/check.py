import pandas as pd
from pandas.core.frame import DataFrame
class Check:
    '''
    result: dict, keys: TP, FP, FN; values: DataFrames
    length: dict, keys: TP, FP, FN; values: ints
    '''
    def __init__(self, df_gold, df_test) -> None:
        self.df_gold = df_gold
        self.df_test = df_test
        self.df_test.index = range(len(df_test))
        self.result = {}
        self.length = {}
    def _zipList(self, df) -> list:
        key = 'mutation' if 'mutation' in df.columns.tolist() else 'processed'
        return list(zip(
            *[list(df['pmid'].astype(str)),list(df[key].replace("X","*", regex=True))]
            ))
    def _getLength(self):
        for k, v in self.result.items():
            self.length[k] = len(v)
class CheckList(Check):
    '''quicker, less info'''
    def __init__(self, df_gold, df_test) -> None:
        super().__init__(df_gold, df_test)
    def compare(self):
        set_gold = set(self._zipList(self.df_gold))
        set_test = set(self._zipList(self.df_test))
        result = {
            'TP': set_gold & set_test,
            'FN': set_gold - set_test,
            'FP': set_test - set_gold
        }
        for k, v in result.items():
            df = pd.DataFrame(list(v), columns = ['pmid', 'mutation'])
            self.result[k] = df.sort_values(by='pmid')
        self._getLength()
class CheckDataFrame(Check):
    '''slower, more info'''
    def __init__(self, df_gold, df_test) -> None:
        super().__init__(df_gold, df_test)
    def compare(self):
        FN = pd.DataFrame()
        gold = self._zipList(self.df_gold)
        test = self._zipList(self.df_test)
        for i in range(len(test)):
            if test[i] in gold: self.df_test.loc[i,"testing"] = "TP"
            else: self.df_test.loc[i,"testing"] = "FP"
        for i in range(len(gold)):
            if gold[i] not in test: 
                FN = FN.append(self.df_gold.loc[i])
            else: continue
        assert isinstance(self.df_test, DataFrame)
        for k in ['TP','FP']:
            new_df = self.df_test[self.df_test.testing == k].copy()
            new_df.drop_duplicates(['pmid', 'processed'], inplace=True)
            new_df.index = range(len(new_df))
            self.result[k] = new_df
        self.result['FN'] = FN
        self._getLength()
    