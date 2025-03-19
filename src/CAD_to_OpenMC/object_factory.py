# SPDX-FileCopyrightText: 2025 E. B. Knudsen <erik@united-neux.eu>
#
# SPDX-License-Identifier: MIT

class ObjectFactory:
    def __init__(self):
        self._builders = {}

    def register_builder(self, key, builder):
        self._builders[key] = builder

    def create(self, key, **kwargs):
        builder = self._builders.get(key)
        if not builder:
            raise ValueError(key)
        return builder(**kwargs)
