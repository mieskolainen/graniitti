# Function call caching mechanism
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


from functools import wraps


def icecache(function):
    """
    Decorate your functions with @icecache
    """
    cache = {}
    @wraps(function)
    def wrapper(*args):
        try:
            return cache[args]
        except KeyError:
            cache[args] = function(*args)
            return cache[args]
    return wrapper

