from abc import ABC, abstractmethod

# TODO: determine which columns are required


class BaseConverter(ABC):
    def __init__(self):
        super().__init__()

    @property
    @abstractmethod
    def converted_data(self):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def load_files(cls, *args):
        return NotImplementedError

    @classmethod
    @abstractmethod
    def convert(cls):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def write(cls, *args):
        raise NotImplementedError

    # NOTE: add functions for additional formats here if we need them

    @staticmethod
    def assert_columns(df):
        # TODO: determine which columns are required
        # TODO: how should we apply this assertion?
        required_columns = []
        assert all([col in list(df) for col in required_columns])
