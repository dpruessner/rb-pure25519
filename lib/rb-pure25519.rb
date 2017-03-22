#
# Finite field algebra 
#

require 'prime'
require 'securerandom'




class Integer
  def ferm_is_prime?
    if self.bit_length < 10
      return Prime.first(200).member? self
    end
    Rb25519::FField.rosetta_mod_exp(2, self-1, self) == 1
  end

  def ferm_ndiv(v)
    a = self / v
    b = self - (a * v)
    [a,b]
  end

  def rb25519_clamp
    v = self & 248
    v &= (127 << (31*8))
    v |= ( 64 << (31*8))
    v
  end

  def to_binary_string
    v = self
    ary = []
    while v > 0
      ary << (v & 0xFF)
      v >>= 8
    end
    ary.pack('c*')
  end

end

class String
  #
  # Convert to a binary fixed size; LSB first
  def to_binary
    v = 0
    self.reverse.each_byte do |byte| 
      v = (v << 8) | byte
    end
    v
  end
  def rb25519_clamp
    bytes = self.each_byte.to_a
    bytes[0] &= 248;
    bytes[31] &= 127;
    bytes[31] |= 64;

    return bytes.pack('c*')
  end
end

class Array
  def to_xz
    [ self[0], self[0].f[1] ]
  end
  def to_xy
    self[0] / self[1]
  end
end




class ECInfinity 
  def self.to_xy
    self
  end
  def inspect
    "#<ECInfinity>"
  end
  def to_s; inspect; end
end



module Rb25519
  class FField
    class FFieldValue 
      attr_accessor :val

      def ==(v)
        @val == v.to_i
      end
      def f
        @field
      end

      def to_i
        @val
      end

      def to_f
        @val.to_f
      end

      def initialize(f, v)
        @field = f
        @val = v.to_i

        if @val >= @field.p  || @val < 0
          @val = @val % @field.p
        end
      end

      def +(v)
        FFieldValue.new(@field, @field.add(self, v))
      end

      def -(v)
        FFieldValue.new(@field, @field.sub(self, v))
      end

      def *(v)
        FFieldValue.new(@field, @field.mul(self, v))
      end

      def inv
        FFieldValue.new(@field, @field.inv(self))
      end

      def /(v)
        FFieldValue.new(@field, @field.div(self, v))
      end

      def **(v)
        FFieldValue.new(@field, @field.exp(self, v))
      end

      def sqrt
        @field.sqrt(self).map{|e| FFieldValue.new(@field, e) }
      end
        
      def -@()
        self * -1
      end

      def inspect
        "#<FFieldValue_#{@field.p}: #{@val} >"
      end

      def to_s
        inspect
      end

    end # FFieldValue




    #
    # Class  methods
    #
    def self.rosetta_mod_exp(b, exp, mod)
      exp < 0 and raise ArgumentError, "negative exponent"
      prod = 1
      base = b % mod
      until exp.zero?
        exp.odd? and prod = (prod * base) % mod
        exp >>= 1
        base = (base * base) % mod
      end
      prod
    end

    def self.eea(i,j)
      s,t,u,v = 1,0,0,1
      while (j != 0)
        q,    r    = i / j,   i % j
        unew, vnew = s    ,   t
       
        s = u - (q * s)
        t = v - (q * t)

        i,  j    =  j,    r
        u,  v    =  unew, vnew


      end
      d, m, n = i, u, v

      return [d, m, n]
    end




    attr_accessor :p

    def initialize(size)
      raise RuntimeError.new("Field must be prime") unless size.ferm_is_prime?
      @p = size
    end


    def [](v)
      FFieldValue.new(self, v)
    end

    def add(a,b)
      nv = a.to_i + b.to_i
      nv -= @p if nv >= @p
      nv
    end

    def sub(a,b)
      nv = a.to_i - b.to_i
      nv += @p if nv < 0
      nv
    end

    def mul(a,b)
      # Naive implementation of multiply
      (a.to_i * b.to_i) % @p
    end

    def inv(v)
      #puts "Inversion"
      return self.class.eea(@p, v.to_i)[1]
    end



    ## 
    # rv = (a / b)
    #
    def div(a,b)
      a.to_i * inv(b.to_i)
    end


    ##
    # rv = b^e
    #
    def exp(b,e)
      self.class.rosetta_mod_exp(b.to_i, e.to_i, @p)
    end

    def sqrt(n)
      n = n.to_i
      return nil if exp(n, (@p-1)/2) != 1

      if (@p % 4) == 3
        r = exp(n, (p+1) / 4)
        return [ r, -r % @p ]
      end



      ##
      ## Implement Tonelli-Shanks (from https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm )
      ## 

      # Factor out Q and S
      #
      #
      p1 = @p-1

      q,s = nil,nil

      p1.bit_length.times.each do |_s|

        _q, _res = p1.ferm_ndiv(1 << _s)
        ## puts [_s, _q, _res].inspect

        if _res == 0 and _q.odd?
          q,s = _q, _s
        end
      end

      ## puts "q,s: #{[q,s].inspect}"

      # Find `z` such that Legendre ( z | p ) == -1
      z = nil
      (1..@p).each{|_z| (z = _z; break) if self.exp(_z, (@p-1)/2) > 1 }
      ## puts "Found z: #{z}"

      c = exp(z, q)
      ## puts "Calculated c: #{c}"

      r = nil
      _r = exp(n, (q+1) / 2 )
      t = exp(n, q)
      m = s


      @p.times do

        if (t == 1)
          r = _r
          ## puts "R is #{r}"
          break
        end

        # Modify t and R
        #

        # find i such that 0 < i < M, such that t**2**i == 1 (mod p)
        i = nil
        i = (1..(m-1)).find{|_i| exp(t, (1<<_i)) == 1 }
        ## puts "Found i: #{i}"

        b = exp(c, (1 << (m - i - 1)))
        _r = mul(_r, b)
        t = mul(t, exp(b, 2) )
        c = exp(b, 2)

        m = i

        ## puts({:b => b, :r => _r, :t => t, :c => c}.inspect)

      end
      

      [r, @p - r] if r
    end



    class EC
      attr_reader :field
      def initialize(field, coeffs=nil)
        @coeffs = coeffs
        @field = field
      end

      def on_curve(x,y)
        raise NotImplementedError.new
      end

      def naive_points
        points = [ECInfinity]
        @field.p.times do |x|
          @field.p.times do |y|
            points << [x,y] if on_curve(x,y)
          end
        end
        points
      end

      def point_add(point_a, point_b)
        xa = point_a[0].kind_of?(FFieldValue) ? point_a[0] : @field[point_a[0]]
        xb = point_b[0].kind_of?(FFieldValue) ? point_b[0] : @field[point_b[0]]

        ya = point_a[1].kind_of?(FFieldValue) ? point_a[1] : @field[point_a[1]]
        yb = point_b[1].kind_of?(FFieldValue) ? point_b[1] : @field[point_b[1]]

        if xa == xb and ya == yb
          return double_point(point_a)
        end
        #puts "point_add: #{point_a.inspect}   + #{point_b.inspect}"

        # All the following operations are in F_p (eg, "mod p")
        
        s = (yb - ya) / (xb - xa)
        #puts "Slope: #{s}"

        xc = s**2 - xa - xb
        yc = (ya * -1) + (xa - xc) * s

        [xc, yc]
      end

      def scale_naive(k, point_a)
        point = point_a

        (k-1).times do
          point = point_add(point, point_a)
        end

        point
      end

      def scale_double_add(k, point_a)
        t = point_a

        bits = k.bit_length

        (bits-1).times.to_a.reverse.each do |bit|
          t = point_add( t, t )
          if (k >> bit) & 0x1 == 1
            t = point_add(t, point_a)
          end
        end
        t
      end

    end


    ##
    #
    # Toy Montgomery curve to test
    #
    # Montgomery curves have the form:
    #
    # B * y^2  =  x^3  +  A * x^2  + x
    #
    #
    #
    # Test cases:
    #
    #   E(@a = 3, @b = 1) over p=47
    #
    #   H = [5,8] on curve.
    #
    #
    class MontgomeryEC < EC
      attr_accessor :a, :b

      def initialize(field, a: nil, b:nil)
        super(field)

        @a = @field[ a || 3]
        @b = @field[ b || 1]

        @a24 = (@a + 2) / @field[4]
      end


      def on_curve(x,y)
        x = @field[x] unless x.kind_of? FFieldValue
        y = @field[y] unless y.kind_of? FFieldValue

        (@b * y**2) == (x**3) + x**2 * @a + x
      end


      ##
      #
      # Add points in affine coordinates
      #
      def point_add(point_a, point_b)

        return point_b if point_a == ECInfinity
        return point_a if point_b == ECInfinity

        xa = point_a[0].kind_of?(FFieldValue) ? point_a[0] : @field[point_a[0]]
        xb = point_b[0].kind_of?(FFieldValue) ? point_b[0] : @field[point_b[0]]

        ya = point_a[1].kind_of?(FFieldValue) ? point_a[1] : @field[point_a[1]]
        yb = point_b[1].kind_of?(FFieldValue) ? point_b[1] : @field[point_b[1]]

        #puts "MontgomeryEC#point-add: #{[xa, ya].inspect}  + #{[xb, yb].inspect}"

        if xa == xb and ya == -yb
          return ECInfinity
        end

        if xa == xb and ya == yb
          return double_point(point_a)
        end

        # All the following operations are in F_p (eg, "mod p")

        l = ( yb - ya) / (xb - xa)
        m = ya - l * xa

        xc = @b * l**2 - @a - xa - xb
        yc = (xa * 2 + xb + @a) * (yb - ya) / (xb - xa)   -   ( @b * (yb - ya) ** 3 ) / (xb - xa)**3     - ya
        [xc, yc]
      end


      ##
      # 
      # Doubles a point in affine coordinates
      #
      def double_point(point_a)
        #puts "Double point: #{point_a.inspect}"

        return point_a if point_a == ECInfinity

        xa = point_a[0].kind_of?(FFieldValue) ? point_a[0] : @field[point_a[0]]
        ya = point_a[1].kind_of?(FFieldValue) ? point_a[1] : @field[point_a[1]]


        bb_inv = (@b * 2 * ya).inv

        c1 = (xa**2 * 3   +    @a * xa * 2     + 1 )

        
        xc = @b * c1**2 * bb_inv**2   - @a - xa - xa
        yc = (xa * 2 + xa + @a) * c1 / (@b * ya * 2)     -    @b * c1**3 * bb_inv**3  - ya


        #x3 = b*   (3*x12+2*a*x1+1) **2   /   (2*b*y1)**2  -a-x1-x1

        #y3 = (2*x1+x1+a)  *(3*x12+2*a*x1+1)/(2*b*y1)-b*(3*x12+2*a*x1+1)3/(2*b*y1)3-y1

        [xc, yc]
      end

    end



    class EC25519 < MontgomeryEC
      attr_accessor :gen

      def initialize
        super(FField.new(2**255 - 19), a: 486662, b: 1 )
      end
    end

  end





  ##

  # Extend the Montgomery EC with projective (XZ) coordinate functions
  #
  class FField::MontgomeryEC

    def xz_from_xy(point)
      return [ @field[point[0].to_i], @field[1] ]
    end


    ##
    # Convert XZ to XY coordinates.
    #
    # point must be Array<FastFrac<FFieldValue>> to work in the EC field.
    #
    #
    # Montgomery curves have the form:
    #
    # B * y^2  =  x^3  +  A * x^2  + x
    #
    def xz_to_xy(point)

      x = point[0] / point[1]

      y_sq = (x**3 + x**2 * @a + x) / @b

      return y_sq.sqrt.map{|e| [x, e]}
    end


    ## 
    #
    # Implementing "mdbl-1987-m" described in DJB's EFD database:
    # http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
    #
    def xz_double(pa)
      x1 = pa[0]
      z1 = pa[1]

      c = (x1 - z1)**2
      d = (x1 * 4 * z1)

      x3 = (x1 + z1) ** 2 *  c
      z3 = d*( c + @a24 * d )
      [x3, z3]
    end

    ###
    ### XZ scaling:
    #
    # Check out: Montgomery Scalar Multiplication for Genus 2 Curves
    # p. 3 (Prop 1) + explanation of the ladder in Efficient Elliptic Curve Point Multiplication
    # 

    ##
    #
    # Implement the XZ coordinates from
    # "Speeding the Pollard and elliptic curve methods of factorization" p.261
    #
    # Actually, best description is at:
    # 
    # Cryptographic Algorithms on Reconfigurable Hardware p.301
    #
    # Actually-- I'm not sure where this one came from.  Figuring out XZ
    # projective point adding was a real pain in the ass!
    # 
    #
    # Test points for scaling/point addition and testing the Montgomery ladder.
    #
    #(5 : 8 : 1)
    # --
    # (2, (14 : 44 : 1))
    # (3, (41 : 36 : 1))
    # (4, (34 : 6 : 1))
    # (5, (23 : 37 : 1))
    # (6, (17 : 4 : 1))
    # (7, (43 : 36 : 1))
    # (8, (8 : 17 : 1))
    # (9, (40 : 28 : 1))
    # (10, (7 : 11 : 1))
    # (11, (46 : 1 : 1))
    # (12, (27 : 29 : 1))
    # (13, (20 : 14 : 1))
    # (14, (6 : 1 : 1))
    # (15, (35 : 14 : 1))
    # (16, (36 : 14 : 1))
    # (17, (45 : 7 : 1))
    # (18, (18 : 17 : 1))
    # (19, (39 : 1 : 1))
    # (20, (37 : 29 : 1))
    # (41, (41 : 11 : 1))
    # (42, (14 : 3 : 1))
    # (43, (5 : 39 : 1))
    # (44, (0 : 1 : 0))
    # (45, (5 : 8 : 1))
    #
    def xz_simple_add(pa, pb, x)

      return pb if (pa == ECInfinity)
      return pa if (pb == ECInfinity)
      
      x1 = pa[0]
      z1 = pa[1]
      x2 = pb[0]
      z2 = pb[1]

      x3 =       ( (x1 - z1)*(x2 + z2) + (x1 + z1)*(x2 - z2) )**2
      z3 =  x *  ( (x1 - z1)*(x2 + z2) - (x1 + z1)*(x2 - z2) )**2

      [x3, z3]
    end


    def scale_proj(k, p)
  #    puts "Scaling #{k} times: #{p.inspect}"

      pa = ECInfinity
      pb = p.to_xz

      x = p[0]


      bits = k.bit_length
  #    puts "Bits: #{bits}"

      (1..bits).each do |j|
  #      puts "Aff[a:x] = #{[pa, pa.to_xy ]}"
  #      puts "Aff[b:x] = #{[pb, pb.to_xy ]}"

        if (k >> (bits - j) ) & 1 == 0

  #        puts "[[ bit: 0 ]]; pb = pa + pb; pa = 2*pa"

          pb = xz_simple_add( pa, pb, x )
          pa = xz_double( pa )
        else

  #        puts "[[ bit: 1 ]]; pb = 2*pb; pa = pa + pb"

          pa = xz_simple_add( pa, pb, x )
          pb = xz_double(pb)
        end

  #      puts

      end

  #      puts "--end--"
  #      puts "Aff[a:x] = #{pa[0] / pa[1]}"
  #      puts "Aff[b:x] = #{pb[0] / pb[1]}"

      return ECInfinity if pa[1] == 0
      pa
    end


    ##
    # 
    # List of scaled points of [5,8] on toy curve to test laddering and other
    # REPL-style exploration/testing to get this working right.
    #
    def pts
      [
        [ @field[ 0], @field[ 0] ],    # 0
        [ @field[ 5], @field[ 8] ],    # 1
        [ @field[14], @field[44] ],    # 2
        [ @field[41], @field[36] ],    # 3
        [ @field[34], @field[ 6] ],    # 4
        [ @field[23], @field[37] ],    # 5
        [ @field[17], @field[ 4] ],    # 6
        [ @field[43], @field[36] ],    # 7
        [ @field[ 8], @field[17] ],    # 8
        [ @field[40], @field[28] ],    # 9
      ]
    end

  end







  # Module Methods

  CURVE   = FField::EC25519.new   
  BASE_XZ = [ CURVE.field[9], CURVE.field[1] ]



  def self.string_to_number(val)
    v = 0
    val.reverse.each_byte do |byte| 
      v = (v << 8) | byte
    end
    v
  end

  def self.number_to_string(v)
    ary = []
    while v > 0
      ary << (v & 0xFF)
      v >>= 8
    end
    ary.pack('c*')
  end

  def self.clamp_string(str)
    bytes = str.each_byte.to_a
    bytes[0] &= 248;
    bytes[31] &= 127;
    bytes[31] |= 64;

    return bytes.pack('c*')
  end

  def self.random_secret_str
    rv = SecureRandom.random_bytes(32)
    rv = clamp_string(rv)
    rv
  end

  def self.random_secret_num
    string_to_number(random_secret_str)
  end
    

  #
  #
  def self.public_key_num(secret)
    if String === secret
      secret = string_to_number(secret)
    end

    xz = CURVE.scale_proj( secret, BASE_XZ )
    (xz[0] / xz[1]).to_i
  end


  def self.public_key_str(secret)
    number_to_string( public_key_num(secret) )
  end

  ##
  # Secret is a 'k' in Q = k*P
  #
  # We want to calculate: 
  #
  #   P_shared_secret = (skey + other_skey) * P_base
  #
  # We get there because:
  #
  #   P_other_pkey = (other_skey) * P_base
  #
  #
  # So continuing to scale P_other_pkey by `skey` will get us to
  # P_shared_secret.  The other party is also doing this calculation; the
  # Abelian group property means this operation is commutative.
  #
  # Note that the Points are all X points in XZ projective space.
  #
  #
  def self.shared_secret_num(pkey, skey)
    if String === pkey
      pkey = string_to_number(pkey)
    end
    if String === skey
      skey = string_to_number(skey)
    end

    shared_xz = CURVE.scale_proj( skey, [ CURVE.field[pkey], CURVE.field[1] ] )

    (shared_xz[0] / shared_xz[1]).to_i    # Final projective -> affine inversion
  end
    
  def self.shared_secret_str(pkey, skey)
    number_to_string( shared_secret_num(pkey, skey) )
  end

end
