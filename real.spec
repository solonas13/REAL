Summary: real
Name: real
Version: 0.0.31
Release: 1
Source0: %{name}-%{version}.tar.bz2
License: GPL version 3
Group: biology
BuildRoot: %{_builddir}/%{name}-root

%description
REAL: An efficient REad ALigner for next generation sequencing reads

%prep
%setup -q
%build
./configure --libdir=%{_libdir} --bindir=%{_bindir} --datarootdir=%{_datadir}
make
%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%files
%defattr(-,root,root)
%{_bindir}/real
%{_bindir}/genpat
%{_bindir}/patsort

