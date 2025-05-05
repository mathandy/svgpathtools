from __future__ import annotations
import hashlib
import struct
from typing import Iterable, Union


def hash_numbers(numbers: Iterable[Union[complex, int, float]]) -> int:
    """Platform-agnostic hashing function

    Args:
        numbers (iterable): iterable of complex (or real) floats and/or ints

    Returns:
        the hash as str
    """

    packed_parts = []
    for n in numbers:
        if isinstance(n, complex):
            packed_parts.append(b'C' + struct.pack('>dd', float(n.real), float(n.imag)))
        elif isinstance(n, float):
            packed_parts.append(b'F' + struct.pack('>d', n))
        else:
            packed_parts.append(b'I' + struct.pack('>q', int(n)))  # 64-bit signed int
    packed = b''.join(packed_parts)

    digest = hashlib.sha256(packed).digest()
    return int.from_bytes(digest, byteorder='big')
