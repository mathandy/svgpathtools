from __future__ import annotations
import hashlib
import struct
import decimal
from typing import Iterable, Union


def float_to_ieee754(f: float) -> bytes:
    """Convert float to IEEE 754 binary64"""
    ctx = decimal.Context(prec=17, rounding=decimal.ROUND_HALF_EVEN)
    d = ctx.create_decimal_from_float(float(f))
    return struct.pack('>d', float(d))


def hash_numbers(numbers: Iterable[Union[complex, int, float, bool]]) -> int:
    """Platform-agnostic hashing function

    Args:
        numbers (iterable): iterable of complex (or real) floats and/or ints

    Returns:
        the hash as str
    """

    packed_parts = []
    for n in numbers:
        if isinstance(n, complex):
            packed_parts.append(b'C' + float_to_ieee754(n.real) + float_to_ieee754(n.imag))
        elif isinstance(n, float):
            packed_parts.append(b'F' + float_to_ieee754(n))
        else:
            packed_parts.append(b'I' + struct.pack('>q', int(n)))
    packed = b''.join(packed_parts)

    digest = hashlib.sha256(packed).digest()
    return int.from_bytes(digest, byteorder='big') % (2**64)
