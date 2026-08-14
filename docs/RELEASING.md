# Releasing FluLens

Every release ships two independent artefacts: the **single
HTML file**, which needs no build, and the **desktop app**, which does.

> Every path and command below is relative to the **repository root**, not to
> `docs/` where this file is.

## The short version

```bash
git tag v1.0.0 && git push origin v1.0.0
```

That runs `.github/workflows/release.yml`. It makes the macOS (universal), Windows,
and Linux builds, plus `flulens.html` and `example_run.zip`. It attaches them to a
**draft** release. Review the draft, then publish it.

A separate push to `main` redeploys the browser version to GitHub Pages.

---

## macOS signing — it decides whether users can open the app

An unsigned `.app` from GitHub carries a quarantine flag. On macOS 15
and later, the old right-click → Open bypass is gone. The user sees *"Apple could
not verify FluLens is free of malware"*. The user must then go to **System Settings →
Privacy & Security → Open Anyway**. At that point, a biologist sends you an email
instead of using the tool. So the setup for signing is worth it.

Signing needs a **Developer ID Application** certificate. An *Apple Development*
certificate is a different thing and does not work for distribution. Check
which one you have with:

```bash
security find-identity -v -p codesigning | grep "Developer ID"
```

### One-time setup

1. **Create the certificate.** Go to
   [developer.apple.com/account/resources/certificates](https://developer.apple.com/account/resources/certificates).
   Add a certificate of type **Developer ID Application**. Download it and open
   it, so it goes into your login keychain.

2. **Create an app-specific password** at
   [appleid.apple.com](https://appleid.apple.com) → Sign-In and Security →
   App-Specific Passwords. Notarisation does not accept your Apple ID password.

3. **Find your Team ID.** It is the ten-character code at the top right of the
   developer portal. You can also use:

   ```bash
   security find-identity -v -p codesigning | grep "Developer ID Application"
   ```

   It is the string in parentheses at the end of the identity name.

### Building signed, locally

```bash
export APPLE_SIGNING_IDENTITY="Developer ID Application: Your Name (TEAMID)"
export APPLE_ID="you@example.com"
export APPLE_PASSWORD="abcd-efgh-ijkl-mnop"   # the app-specific password
export APPLE_TEAM_ID="TEAMID"
cd desktop && cargo tauri build --target universal-apple-darwin
```

Tauri signs the `.app`, submits it to notarisation, waits, and staples the ticket.

### Tauri does not notarise the DMG — you must do it

This is the trap. Tauri notarises the **app**. Then it builds the DMG *from* the app
and only **signs** the image. So the app comes out as `Notarized Developer ID`, but
the disk image comes out as `Unnotarized Developer ID, rejected`. A downloaded DMG is
quarantined, and Gatekeeper checks the *image*, not only the app inside it. So if you
skip this step, your users get the exact warning that you signed to prevent. And
nothing in the build output shows that a problem exists.

```bash
DMG=desktop/target/universal-apple-darwin/release/bundle/dmg/FluLens_1.0.0_universal.dmg
xcrun notarytool submit "$DMG" \
  --apple-id "$APPLE_ID" --team-id "$APPLE_TEAM_ID" --password "$APPLE_PASSWORD" --wait
xcrun stapler staple "$DMG"
```

`.github/workflows/release.yml` does this automatically after `tauri-action`.

### Verifying

Check the image and the app separately. A pass on one says nothing about the
other:

```bash
xcrun stapler validate "$DMG"
spctl -a -vvv -t open --context context:primary-signature "$DMG"   # the image
spctl -a -vvv -t install desktop/target/universal-apple-darwin/release/bundle/macos/FluLens.app
```

All must report **accepted** with `source=Notarized Developer ID`.

The correct test is to simulate a download. macOS judges an unquarantined local file
more leniently than a file that came from the internet:

```bash
cp "$DMG" /tmp/dl.dmg
xattr -w com.apple.quarantine "0081;0;Safari;$(uuidgen)" /tmp/dl.dmg
spctl -a -vvv -t open --context context:primary-signature /tmp/dl.dmg
```

### Storing the app-specific password

Keep it in the keychain, not in your shell history or a dotfile:

```bash
security add-generic-password -a "you@example.com" -s FLULENS_NOTARY -w
export APPLE_PASSWORD="$(security find-generic-password -s FLULENS_NOTARY -w)"
```

### Building signed, in CI

Set these repository secrets (Settings → Secrets and variables → Actions):

| Secret | What it is |
|---|---|
| `APPLE_CERTIFICATE` | the Developer ID cert exported as `.p12`, base64-encoded |
| `APPLE_CERTIFICATE_PASSWORD` | the password you set when you exported the `.p12` |
| `APPLE_SIGNING_IDENTITY` | `Developer ID Application: Your Name (TEAMID)` |
| `KEYCHAIN_PASSWORD` | any string; CI makes a temporary keychain with it |
| `APPLE_ID` | your Apple ID email |
| `APPLE_PASSWORD` | the app-specific password |
| `APPLE_TEAM_ID` | the ten-character team ID |

Export and encode the certificate with:

```bash
security find-certificate -c "Developer ID Application" -p   # confirm it exists
# then export from Keychain Access as Certificate.p12, and:
base64 -i Certificate.p12 | pbcopy
```

The workflow runs without these secrets, but it then makes an unsigned build. Know
this in advance, so you do not need to debug it later: a green CI run does **not** mean
the app opens cleanly.

## Windows

The build is unsigned, and SmartScreen warns on first run. An OV/EV
code-signing certificate costs a few hundred dollars a year. Unless Windows users
ask for it, document the warning and do not buy one.

## Version numbers

Three places must agree, and nothing checks them for you:

- `desktop/tauri.conf.json` → `version`
- `desktop/Cargo.toml` → `version`
- `CITATION.cff` → `version` and `date-released`

The git tag must match. `tauri-action` takes the release name from the tag,
not from the config. So a mismatch appears as a release whose title and whose
installer disagree.

## After publishing

Zenodo makes a DOI for each GitHub release, after you enable the repository at
[zenodo.org/account/settings/github](https://zenodo.org/account/settings/github).
Add the DOI to `CITATION.cff` and to the README badge. Flumina already works
this way.
