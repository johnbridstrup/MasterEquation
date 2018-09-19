


class TimeSeriesError(Exception):
    pass


class TimeSeries:
    def __init__(self, data, t_data=None, y_unit=None, t_unit='s'):
        if t_data is None:
            self.t_data = data.pop('t')
            self.data = data
        else:
            self.t_data = t_data
            self.data = data

        self.y_unit = y_unit
        self.t_unit = t_unit


