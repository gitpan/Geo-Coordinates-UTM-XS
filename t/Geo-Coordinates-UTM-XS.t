#!/usr/bin/perl

use Test::More tests => 9991;
BEGIN { use_ok('Geo::Coordinates::UTM::XS') };

use constant maxerror => 1e-2;

use warnings;
use strict;

sub fleq ($$;$) {
    if (abs($_[0] - $_[1]) < maxerror) {
        pass($_[2]);
    }
    else {
        fail($_[2]);
        diag("floating point value $_[0] too different to reference $_[1]");
    }
}

open TD, '<', 't/test.dat'
    or die "unable to open test data file 't/test.dat'";

my $latlon = "CCDEFGHJKLMNPQRSTUVWXX";

while(<TD>) {
    chomp;
    next if /^\s*(?:#.*)?$/;
    my ($ellipsoid, $latitude, $longitude, $zone, $easting, $northing) = split /\|/;
    my ($z, $e, $n) = latlon_to_utm($ellipsoid, $latitude, $longitude);
    is($z, $zone, "zone $.");
    fleq($e, $easting, "easting $.");
    fleq($n, $northing, "northing $.");

    my ($lat, $lon) = utm_to_latlon($ellipsoid, $z, $easting, $northing);
    fleq($lon, $longitude, "longitude $.");
    fleq($lat, $latitude, "latitude $.");

    my ($zone_number, $zone_letter) = $zone =~ /^(\d+)(\w)/;
    ($z, $e, $n) = latlon_to_utm_force_zone($ellipsoid, $zone_number, $latitude, $longitude);
    is($z, $zone, "fz zone $.");
    fleq($e, $easting, "fz easting $.");
    fleq($n, $northing, "fz northing $.");



    my $z1 = $zone_number + int(-2 + rand 5);
    $z1 -= 60 if $z1 > 60;
    $z1 += 60 if $z1 < 1;

    my $l1 = ($latlon =~ /(.)($zone_letter)(.)/, '')[rand(4)];
    ($z, $e, $n) = latlon_to_utm_force_zone($ellipsoid, "$z1$l1", $latitude, $longitude);
    ($lat, $lon) = utm_to_latlon($ellipsoid, $z, $e, $n);
    fleq($lon, $longitude, "fz longitude (zone $zone) $.");
    fleq($lat, $latitude, "fz latitude (zone $zone) $.");
}

