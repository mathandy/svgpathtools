import hashlib
import struct
import sys


def hash_numbers(numbers):
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
    return int.from_bytes(digest, byteorder='big') if sys.version_info[0] >= 3 else long(int(digest.encode('hex'), 16))
