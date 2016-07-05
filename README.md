# Ruby Curve25519

DJB's [Curve25519](http://cr.yp.to/ecdh.html)  implemented in pure Ruby using
core Bignum for calculations on prime Finite Fields.

## Usage

### Ephermeral DH Exchange

```ruby
require 'rb-pure25519'
require 'securerandom'

secret = SecureRandom.random_bytes(32).rb25519_clamp

pkey = Rb25519
```



Powered by Curve25519.
