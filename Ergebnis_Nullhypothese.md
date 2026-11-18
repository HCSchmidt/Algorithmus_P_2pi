# Auswertung der P(2π)-Massenformel – Ergebnisse und Fazit

## Kurzfassung

Mit den Daten in diesem Repo gibt es keinen statistischen Beleg, dass P(2π) Teilchenmassen beschreibt. Was die Formel trifft, erklärt sich durch zwei Dinge: ein sehr dichtes Spektrum und nachträglich angepasste Konstanten.

## 1. Der Code im Repo rechnet nicht die Formel aus dem Paper

Ich habe die Python-Funktion aus dem Paper (Abb. 3) abgeschrieben und mit `energy_model.py` verglichen. Bei 99,85 % aller Koeffizienten-Kombinationen kommen unterschiedliche Werte heraus, mit Abweichungen bis 0,26 mₑ.

- **Ursache:** der Term E7. Im Paper wird er innerhalb der Paarungsschleife abgezogen, und zwar nur für das erste Null/Null-Paar (danach `break`). Das Repo zieht ihn danach für alle Null/Null-Paare ab.
- **Folge:** Mit dem Repo-Code ergibt das Proton-Tupel 1835,894 statt 1836,1527. Alle Ergebnisse in `results/` stammen also von einer anderen Formel als der im Paper.

Mit der Paper-Formel kommen die Werte aus dem Paper exakt heraus (Proton, Neutron, Myon, Pionen).

## 2. Fast alle „Treffer“ sind Zufall

Das ist der Kern der Auswertung. Alle Kombinationen zusammen ergeben etwa 4 Millionen positive Energiewerte zwischen 0 und 3718 mₑ, also rund 1000 Werte pro mₑ. Jedes Messfenster, das breiter als etwa 0,001 mₑ ist, wird deshalb zwangsläufig getroffen, egal welche Masse man hineinlegt.

Für jedes Teilchen habe ich berechnet, wie viele Treffer man allein durch Zufall erwarten würde: λ = Dichte × Fensterbreite.

- **Ergebnis:** Bei 31 getroffenen Teilchenfenstern stimmt die Zahl der Treffer praktisch genau mit λ überein (Median-Verhältnis 1,00). Beispiel s-Quark: 30.101 Treffer, erwartet 28.680.
- **Gegenprobe mit erfundenen Massen:** Ich habe 20.000 zufällige Massen mit gleich breiten Fenstern getestet. Davon werden 64,9 % getroffen, der Zufall sagt 65,3 % voraus.

Die „Wahrscheinlichkeiten“ Δi/2π im Paper (z. B. s = 49 %) messen deshalb nur, wie dicht das Spektrum ist. Über Physik sagen sie nichts.

## 3. Nur drei Treffer schlagen den Zufall: Proton, Neutron, Myon

Genau diese drei Teilchen wurden zur Konstruktion der Konstanten benutzt:

- E_C₊ (8 Terme mit π⁻ᵏ) für das Proton
- E₀ = π⁻¹² + 2π⁻¹⁴ für das Neutron
- E_C₋ = 2π − π⁻¹ + E_C₊ für das Myon

**Test der Anpassungsfreiheit:** Ich nehme einen Gitterwert und erlaube eine Korrektur der gleichen Art, also Σ aₖπ⁻ᵏ mit k = 7…14 und aₖ ∈ {−2…2}. Damit treffe ich 93,5 % beliebiger erfundener Zielmassen mit Proton-Genauigkeit, also 10 Stellen. Mit anderen Basen statt 2π sind es immerhin 69–87 %.

Die 10-stellige Übereinstimmung beim Proton zeigt also vor allem, wie gut die Konstante angepasst wurde. Ein Beleg für die Theorie ist sie nicht. Details dazu in Abschnitt 9.

## 4. Andere Basen statt 2π

Ich habe π überall durch b/2 ersetzt und 206 andere Basen getestet, darunter 6, 7, 2e, √40 und 200 zufällige Werte. Keine davon trifft Proton oder Neutron.

Das klingt zunächst beeindruckend, ist aber nur logisch: Die Konstanten wurden für 2π gebaut. Wie Abschnitt 3 zeigt, könnte man sie für jede andere Basis genauso bauen.

## 5. Der eigentliche Test: Teilchen, die nicht eingebaut wurden

Aussagekräftig sind nur Massen, die sehr genau gemessen sind und nicht zur Konstruktion verwendet wurden. Alle werden verfehlt:

| Teilchen | Ergebnis |
|---|---|
| Elektron (genau 1 mₑ per Definition) | verfehlt, beste Näherung 0,9988 |
| Wasserstoffatom ¹H | verfehlt |
| Myonium | verfehlt |
| Positronium | verfehlt |
| Deuteron | verfehlt |
| Deuterium-Atom | verfehlt |

Eine echte Theorie müsste genau hier treffen.

## 6. Die eine interessante Beobachtung

Das „nackte“ Polynom

(2π)⁴ + (2π)³ + (2π)² − 2π − 1 − (2π)⁻¹ = 1838,6317

liegt nur 2,8·10⁻⁵ relativ neben der Neutronenmasse. Für ein festes Ziel ist das ein etwa 0,6 % wahrscheinlicher Zufall.

Das ist vergleichbar mit Lenz' bekannter Beobachtung von 1951: 6π⁵ ≈ mₚ/mₑ. Rechnet man ein, wie viele Formeln und Ziele man hätte ausprobieren können, ist es nicht signifikant. Trotzdem ist es der einzige Befund, der ein „hm, nett“ verdient.

## 7. Weitere Fehler im Code (nicht behoben)

- **Sektor-Schnitte wirken nie:** Der Code vergleicht mit `"heavy"` und `"nucleon"`, die Sektoren heißen aber `"222-c_tau"` und `"H-Atom"`.
- **Bin-Index überlappt:** `256+32·i4+4·i3+i2` wirft ab Sektor 011 verschiedene Polynome in denselben Bin.
- **Wasserstoff-Masse falsch:** 1837,47 ist das Atomgewicht von natürlichem Wasserstoff. Reiner Wasserstoff (¹H) hat 1837,153.
- **Nicht erreichbar:** Ds⁺, Higgs und Top-Quark liegen außerhalb des Gitters (Maximum 3718 mₑ). „Alle Elementarteilchen“ stimmt also nicht.
- **`fuer_papa.py`** ignoriert die eigene `ARGS`-Einstellung.
- **Fehler im Paper selbst:**
  - Die Tupel-Bezeichnungen für Neutron und Proton sind vertauscht.
  - Ein Myon-Tupel ist falsch.
  - Das Vorzeichen von 2π⁻¹⁴ ist uneinheitlich.
  - In der Zusammenfassung steht die rekursive Definition „E_C₋ = 2π − π⁻¹ + E_C₋“.

## 8. Was am Code geändert wurde

- **Schnellere Berechnung:** Beim schweren Sektor dauert die Rechnung jetzt 5 s statt etwa 75 s. Alle 20 Ausgabedateien (CSV, TXT, PNG) sind byte-identisch zu vorher, abgesichert durch Tests.
- **Neue Tests:** 20 Tests laufen durch. 3 sind bewusst als „erwarteter Fehlschlag“ markiert. Sie dokumentieren die Abweichung von Repo und Paper und schlagen erst an, wenn jemand sie behebt.
- **Neue Auswertung:** `make eval` startet sie, den ausführlichen Bericht findest du in `results/evaluation/REPORT.md`.

## 9. Was heißt „Konstanten, die genau für diese Teilchen gebaut wurden“?

Jede Energie in der Formel besteht aus zwei Teilen:

**E = Polynom-Teil (aus den Koeffizienten i₄…i₋₁) + Konstante (abhängig von C)**

Der Polynom-Teil ist das eigentliche „2π-Gitter“. Die Konstante ist ein fester Zusatzwert, der nicht aus einer Regel folgt, sondern im Paper von Hand festgelegt wurde. Bei allen drei Präzisionstreffern entscheidet genau diese Konstante, ob getroffen wird.

### Neutron: E₀ = π⁻¹² + 2π⁻¹⁴

Das Neutron hat das Tupel (1,1,1,1,1,1) mit C = 0. Bei C = 0 wird der Standardwert E₀ = π⁻¹² + 2π⁻¹⁴ ≈ 0,0000013 addiert.

| | Wert | Abweichung vom Messwert |
|---|---|---|
| Messwert | 1838,683662 ± 0,00000074 | |
| Polynom ohne E₀ | 1838,6836611 | −0,00000085 → **knapp daneben** |
| Polynom mit E₀ | 1838,6836624 | +0,00000045 → **Treffer** |

Die Konstante ist winzig, liegt aber genau in der Größenordnung der Messunsicherheit des Neutrons. Ohne sie wäre es kein Treffer. Warum ausgerechnet π⁻¹² + 2π⁻¹⁴ und nicht etwa π⁻¹³? Dafür gibt es keine Regel, außer dass es zum Neutron passt.

### Proton: E_C₊ ist eine schrittweise Annäherung

Das Proton hat dasselbe Tupel wie das Neutron, aber C = 1. Dadurch wird E₀ durch E_C₊ ersetzt. Damit das Proton herauskommt, muss E_C₊ genau −2,5309877 sein. Das ist der Massenunterschied zwischen Neutron und Proton.

So ist E_C₊ aufgebaut, Term für Term:

| Term | Summe | Rest zum Ziel |
|---|---|---|
| −π | −3,14159 | +0,61 |
| +2π⁻¹ | −2,50497 | −0,026 |
| −π⁻³ | −2,53722 | +0,0062 |
| +2π⁻⁵ | −2,53069 | −0,0003 |
| −π⁻⁷ | −2,53102 | +0,00003 |
| +π⁻⁹ | −2,5309864 | −0,0000013 |
| −π⁻¹² | −2,5309875 | −0,0000002 |
| −2π⁻¹⁴ | −2,53098774 | **+0,00000002** ✓ |

Man sieht das Muster: Jeder Term ist kleiner als der vorige, und sein Vorzeichen wechselt immer so, dass der verbleibende Rest kleiner wird. Das funktioniert wie das Aufschreiben einer Zahl in Dezimalstellen, nur zur Basis π statt 10.

Mit diesem Verfahren kann man **jede beliebige Zahl** auf beliebig viele Stellen darstellen. Hätte das Proton 1836,4 gewogen, hätte man einfach andere Vorzeichen und Exponenten gewählt. Dass das Proton auf 10 Stellen stimmt, beweist deshalb nichts: Die 10 Stellen wurden über die 8 Terme hineingesteckt.

Der Test aus Abschnitt 3 zeigt das: Mit einer Konstanten dieser Bauart treffe ich 93,5 % beliebiger, frei erfundener Massen auf Proton-Genauigkeit.

### Myon: E_C₋ = 2π − π⁻¹ + E_C₊

Das Myon hat C = −1 und bekommt damit eine dritte Konstante E_C₋. Sie wird aus E_C₊ plus zwei weiteren Termen gebildet:

- Polynom-Teil: 203,3344
- benötigt bis zum Messwert: +3,4338872
- E_C₋ liefert: +3,4338877 → Abweichung 0,0000005 bei einer Messunsicherheit von 0,0000046 → Treffer

Hier wurden also zwei weitere Bausteine (2π und −π⁻¹) gewählt, damit die Summe passt.

### Warum das entscheidend ist

Eine Theorie ist überzeugend, wenn sie etwas **vorhersagt**, das man vorher nicht hineingesteckt hat. Hier läuft es umgekehrt: Man kennt die Messwerte von Proton, Neutron und Myon und baut die Konstanten so, dass sie herauskommen. Das ist eine **Anpassung**, keine Vorhersage. Physiker sprechen von „Fitten“.

Daher kommt der Test mit den nicht eingebauten Teilchen (Abschnitt 5). Auch das Wasserstoffatom, Positronium, das Deuteron oder das Elektron (genau 1) sind extrem genau gemessen. Für sie wurde keine eigene Konstante gebaut, und alle werden verfehlt.

Das ist das typische Bild einer angepassten Formel: Sie stimmt dort perfekt, wo sie angepasst wurde, und trifft sonst nur im Rahmen des Zufalls.

## 10. Wie wahrscheinlich ist es, dass doch etwas dran ist?

Kurz gesagt: **sehr unwahrscheinlich.** Dass diese Formel einen echten physikalischen Zusammenhang beschreibt, halte ich für deutlich unter 1 %, eher im Promille-Bereich. Das ist eine Einschätzung, keine berechnete Zahl.

### Was dagegen spricht

**1. Die Daten liefern keinen Beleg.** Das ist der wichtigste Punkt. Eine Theorie gewinnt an Glaubwürdigkeit, wenn sie etwas trifft, das sie ohne echten Zusammenhang nicht treffen würde. Hier passiert genau das Gegenteil:

- Die breiten Treffer entsprechen exakt dem Zufall.
- Die drei Präzisionstreffer wurden über die Konstanten eingebaut.
- Jede unabhängige Präzisionsprüfung geht daneben.

Die Auswertung verschiebt die Wahrscheinlichkeit also nicht nach oben, eher etwas nach unten.

**2. Die Physik passt nicht dazu.** Im Standardmodell sind viele dieser Massen keine festen Naturzahlen:

- **Quarkmassen** hängen von der Energieskala und vom Rechenverfahren ab. Die „u-Quark-Masse 2,16 MeV“ gilt bei 2 GeV im sogenannten MS-bar-Schema, bei einer anderen Skala ist sie eine andere Zahl. Eine Formel, die *eine* feste Zahl dafür liefert, beschreibt also eine Konvention, keine Naturkonstante.
- **Die Protonmasse** entsteht zu etwa 99 % aus der Bindungsenergie der starken Wechselwirkung, nicht aus den Quarkmassen. Sie lässt sich heute auf Supercomputern (Gitter-QCD) auf etwa 1 % genau aus der Theorie berechnen. Eine einfache π-Formel müsste diese ganze Dynamik zufällig nachbilden.

**3. Die Vorgeschichte.** Es gibt eine lange Reihe solcher Zahlenspiele:

- Eddington behauptete, die Feinstrukturkonstante sei genau 1/136, später 1/137.
- Lenz bemerkte 6π⁵ ≈ mₚ/mₑ.
- Dazu kommen unzählige ähnliche Formeln aus π, e und √2.

Keine davon hat sich als echter Zusammenhang herausgestellt. Bei so vielen Versuchen entstehen zwangsläufig verblüffende Übereinstimmungen. Die einzige, die bis heute ernsthaft diskutiert wird, ist die Koide-Formel für die Massen der drei geladenen Leptonen, und auch sie ist nicht anerkannt.

**4. Die Formel hat zu viele Freiheiten.** 7 Koeffizienten, 5 C-Werte, Paarungsregeln, handgewählte Konstanten und Sektoren nach Wahl: Mit so vielen Stellschrauben kann man fast alles darstellen. Je mehr eine Formel darstellen *kann*, desto weniger sagt es aus, wenn sie etwas darstellt.

### Was ich fair einräume

- **Die Neutron-Beobachtung ist hübsch.** Das schlichte Polynom (2π)⁴ + (2π)³ + (2π)² − 2π − 1 − (2π)⁻¹ landet auf 2,8·10⁻⁵ relativ an der Neutronenmasse. Einzeln betrachtet ist das ein 0,6 %-Zufall. Es bleibt aber eine einzige Übereinstimmung, ohne Vorhersage und ohne Erklärung, warum gerade diese Kombination.
- **Die Auswertung prüft nur, was im Repo und im Paper steht.** Falls es ungeschriebene Regeln gibt, nach denen C und die Korrekturterme *zwingend* folgen, statt gewählt zu werden, wäre das ein neuer Ausgangspunkt. Davon habe ich aber nichts gefunden.

### Was mich umstimmen würde

Ich würde meine Einschätzung deutlich nach oben korrigieren, wenn alle drei Bedingungen erfüllt sind:

1. **Feste Regeln:** Es gibt eine Vorschrift, die *vorher* festlegt, welches Tupel und welche Konstante zu welchem Teilchen gehört, ohne Auswahl nach dem Ergebnis.
2. **Echte Vorhersagen:** Mit genau diesen Regeln trifft die Formel mehrere hochpräzise Werte, die nicht hineingesteckt wurden, zum Beispiel ¹H, Positronium, Myonium oder mₙ − mₚ, jeweils innerhalb der Messunsicherheit.
3. **Bestätigung durch neue Messungen:** Am stärksten wäre es, wenn die Formel einen Wert vorhersagt, *bevor* er genauer gemessen wird, und die neue Messung ihn bestätigt.

Das sind auch genau die Tests, die sich anbieten, wenn die Idee weiterverfolgt werden soll. Die Werkzeuge dafür liegen im Repo (`make eval`), und neue Kandidaten lassen sich in `src/polynome2pi/evaluation/nulltest.py` als zusätzliche „holdout“-Teilchen eintragen.

## Fazit

Die Formel erzeugt ein so dichtes Zahlenspektrum, dass sie fast jede Masse „trifft“. Die drei hochpräzisen Treffer gehen auf Konstanten zurück, die genau für diese Teilchen gebaut wurden. Bei unabhängigen Präzisionstests versagt sie ausnahmslos.

Damit die Theorie ernst genommen werden kann, müsste sie vorher festgelegte Regeln für C und die Korrekturterme haben. Mit diesen Regeln müsste sie dann Massen vorhersagen, die nicht in ihre Konstruktion eingeflossen sind, etwa Wasserstoff, Positronium oder mₙ − mₚ.

Die Wahrscheinlichkeit, dass hier dennoch ein echter physikalischer Zusammenhang vorliegt, schätze ich auf deutlich unter 1 % (Abschnitt 10).
