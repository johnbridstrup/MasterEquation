import pandas as pd
import holoviews as hv
from bokeh.sampledata import stocks
from holoviews.operation.timeseries import rolling, rolling_outlier_std
from holoviews.streams import Stream
hv.notebook_extension('bokeh')

def load_symbol(symbol, **kwargs):
    df = pd.DataFrame(getattr(stocks, symbol))
    df['date'] = df.date.astype('datetime64[ns]')
    return hv.Curve(df, ('date', 'Date'), ('adj_close', 'Adjusted Close'))

stock_symbols = ['AAPL', 'FB', 'GOOG', 'IBM', 'MSFT']
dmap = hv.DynamicMap(load_symbol, kdims='Symbol').redim.values(Symbol=stock_symbols)